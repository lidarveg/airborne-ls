#!/usr/bin/env python

"""
las/laz file manipulations for lidarveg processing of airborne lidar

key activites / options:

-- check file for issues with supplied file
-- build a processing index
-- standardise laz version converting from las where required
-- update projections where needed

"""

import logging
import os
import struct
import sys
import zipfile
from copy import copy
from pathlib import Path

import laspy
import numpy as np

from airborne_ls import filenaming_methods, qvf, const

logger = logging.getLogger(__name__)


class LazNdx:
    """
    Hold the data from our simple square-bin index, as read in from the VLRs
    """
    def __init__(self, filehandle):
        """
        Read the index from the open file

        Parameters:
          filehandle (laspy.lasreader.LasReader): Open LAS file
        """
        binInfoVlr = binBoundsVlr = None

        vlrs = filehandle.header.vlrs
        binInfoList = vlrs.get_by_id(user_id=const.VLR_USERID_JRSRP,
                                     record_ids=[const.VLR_RECORDID_BINSIZE])
        if len(binInfoList) > 0:
            binInfoVlr = binInfoList[0]
            binBoundsList = vlrs.get_by_id(user_id=const.VLR_USERID_JRSRP,
                                           record_ids=[const.VLR_RECORDID_BINBOUNDS])
            binBoundsVlr = binBoundsList[0]
            oldVlrs = False
        else:
            # Old non-conformant VLRs
            for vlr in vlrs:
                if vlr.user_id == "BINSIZE":
                    # VLR containing binSize and nbins
                    binInfoVlr = vlr
                elif vlr.user_id == "BIN_POS":
                    # VLR containing bin positions
                    binBoundsVlr = vlr
            oldVlrs = True

        if binInfoVlr is None or binBoundsVlr is None:
            raise ValueError("Unable to find bin index VLRs in file")

        # Unpack binSize and nbins
        binInfoFormat = "<dQ"
        if oldVlrs:
            binInfoFormat = "2d"
        binSize, nBins = struct.unpack(binInfoFormat, binInfoVlr.record_data)

        if nBins > 1:
            numElements = int(nBins + 1)
            binIndexFormat = f"<{numElements}Q"
            if oldVlrs:
                binIndexFormat = f"{numElements}d"
            binBounds = struct.unpack(binIndexFormat, binBoundsVlr.record_data)
            binBounds = np.array(binBounds, dtype=np.int32)
        else:
            raise ValueError("Failed to import laz binning index")

        self.binSize = binSize
        self.nBins = nBins
        self.binBounds = binBounds


class BinnedData:
    """
    Hold the point data from a LAZ file, divided into bins as-per the index
    """
    def __init__(self, filename, row=None, col=None):
        """
        Open the file, and read data into binned structure. Limit to nominated row or
        col if required. If both row and col are given, then only the intersecting bin
        is loaded.
        """
        f = laspy.open(filename)
        self.ndx = LazNdx(f)
        self.dataByBin = {}
        # We assume that tiles are always square, and thus same number of bin rows as cols
        self.numBinRows = int(np.round(np.sqrt(self.ndx.nBins)))
        self.numBinCols = self.numBinRows

        # We read the whole dataset in. I would love to be able to just read in the points
        # we need, but the compression scheme makes this complicated. Doing it the obvious
        # way ends up taking longer than just reading the whole lot in.
        # The data is in the form of our numpy recarray structure.
        data = self.laspy2rec(f.read())

        # Split the data into the desired bins
        if row is None and col is None:
            for r in range(self.numBinRows):
                for c in range(self.numBinCols):
                    binNum = self.makeBinNum(r, c)
                    i1 = self.ndx.binBounds[binNum]
                    i2 = self.ndx.binBounds[binNum + 1]
                    self.setData(r, c, data[i1:i2])
        elif row is not None and col is None:
            for c in range(self.numBinCols):
                binNum = self.makeBinNum(row, c)
                i1 = self.ndx.binBounds[binNum]
                i2 = self.ndx.binBounds[binNum + 1]
                self.setData(row, c, data[i1:i2])
        elif col is not None and row is None:
            for r in range(self.numBinRows):
                binNum = self.makeBinNum(r, col)
                i1 = self.ndx.binBounds[binNum]
                i2 = self.ndx.binBounds[binNum + 1]
                self.setData(r, col, data[i1:i2])
        else:
            binNum = self.makeBinNum(row, col)
            i1 = self.ndx.binBounds[binNum]
            i2 = self.ndx.binBounds[binNum + 1]
            self.setData(row, col, data[i1:i2])

        del data

    @staticmethod
    def laspy2rec(data):
        """
        Convert from laspy's point record laspy.lasdata.LasData strusture into our own
        homegrown numpy recarray. This holds only the columns we want to use, and is
        more memory-efficient. It also allows copying of sub-arrays, allowing us to
        easily split into bins.

        Parameters:
          data (LasData): Points read directly from file

        Returns:
          (numpy recarray): Custom record array of point data
        """
        recarray = np.rec.fromarrays(
            [
                data.return_num,
                data.num_returns,
                data.gps_time,
                data.intensity,
                data.classification,
                data.x,
                data.y,
                data.z,
            ],
            names=[
                "RETURN_NUMBER",
                "NUMBER_OF_RETURNS",
                "TIMESTAMP",
                "INTENSITY",
                "CLASSIFICATION",
                "X",
                "Y",
                "Z",
            ],
            formats=["u1", "u1", "<f8", "<i4", "u1", "<f8", "<f8", "<f8"],
        )
        return recarray

    @staticmethod
    def makeBinName(r, c):
        """
        Return internal bin name for given bin row/col.
        Top-left bin is (r, c) == (0, 0)
        """
        return f"row_{r}_col_{c}"

    def makeBinNum(self, r, c):
        """
        Return internal bin number for given bin row/col.
        Top-left bin is (r, c) == (0, 0)
        """
        return r * self.numBinCols + c

    def getData(self, row, col):
        """
        Return the point data for the nominated bin.
        Top-left bin is (r, c) == (0, 0)
        """
        binName = self.makeBinName(row, col)
        return self.dataByBin.get(binName, None)

    def setData(self, r, c, data):
        """
        Set the given point data for the nominated bin row/col.
        Top-left bin is (r, c) == (0, 0)

        Parameters:
          r, c (int): Row/col of nominated bin (row zero is top row)
          data (numpy recarray): Point data for this bin
        """
        binName = self.makeBinName(r, c)
        self.dataByBin[binName] = np.copy(data)

    def numBins(self):
        """
        Return total number of bins stored in this instance
        """
        return len(self.dataByBin)

    def numPoints(self):
        """
        Return total number of points stored in this instance (all points for
        all bins)
        """
        return sum([len(self.dataByBin[k]) for k in self.dataByBin])

    def mergeNeighbour(self, other, direction):
        """
        Merge the given 'other' instance into the current one. Direction is from the perspective
        of the central tile, so 'NW' means that 'other' lies to the north-west of the central tile
        We copy over just the immediate next bin(s), according to direction.
        """
        numBinRows = self.numBinRows
        numBinCols = self.numBinCols
        otherNrows = other.numBinRows
        otherNcols = other.numBinCols
        if otherNrows != numBinRows or otherNcols != numBinCols:
            raise ValueError("Bin row/col count mis-match")

        if direction == 'SW':
            self.setData(numBinRows, -1, other.getData(0, otherNcols - 1))
        elif direction == 'S':
            for c in range(numBinCols):
                self.setData(numBinRows, c, other.getData(0, c))
        elif direction == 'SE':
            self.setData(numBinRows, numBinCols, other.getData(0, 0))
        elif direction == 'W':
            for r in range(numBinRows):
                self.setData(r, -1, other.getData(r, otherNcols - 1))
        elif direction == 'E':
            for r in range(numBinRows):
                self.setData(r, numBinCols, other.getData(r, 0))
        elif direction == 'NW':
            self.setData(-1, -1, other.getData(otherNrows - 1, otherNcols - 1))
        elif direction == 'N':
            for c in range(numBinCols):
                self.setData(-1, c, other.getData(otherNrows - 1, c))
        elif direction == 'NE':
            self.setData(-1, numBinCols, other.getData(otherNrows - 1, 0))


def readBinnedData(filename, tileSize=None, withNeighbours=False):
    """
    Read point data from the given LAZ file, divide into bins according to the stored index

    Parameters:
      filename (str): Name of LAZ file to read
      tileSize (float): Size (i.e. length of side) (metres) of tile in the file
      withNeighbours (bool): If True, also read data from neighbouring bins in
                             all 8 neighbouring tiles (wherever available). Assumes
                             our standard internal file naming.

    Returns:
      (BinnedData): The data for the tile, binned by the index
    """
    if withNeighbours and (tileSize is None):
        raise ValueError("withNeighbours requires tileSize")

    binnedData = BinnedData(filename)
    if withNeighbours:
        offsetByDirection = {
            'SW': (-tileSize, -tileSize), 'S': (0, -tileSize), 'SE': (tileSize, -tileSize),
            'W': (-tileSize, 0), 'E': (tileSize, 0),
            'NW': (-tileSize, tileSize), 'N': (0, tileSize), 'NE': (tileSize, tileSize)
        }
        rowColToRequest = {
            'SW': (0, binnedData.numBinCols - 1), 'S': (0, None), 'SE': (0, 0),
            'W': (None, binnedData.numBinCols - 1), 'E': (None, 0),
            'NW': (binnedData.numBinRows - 1, binnedData.numBinCols - 1),
            'N': (binnedData.numBinRows - 1, None), 'NE': (binnedData.numBinRows - 1, 0)
        }
        where = qvf.getwhere(filename)
        for direction in offsetByDirection:
            (xOffset, yOffset) = offsetByDirection[direction]

            nbrWhere = filenaming_methods.neighbourTileWhere(where, xOffset, yOffset)
            nbrFilename = qvf.setwhere(filename, nbrWhere)
            if os.path.exists(nbrFilename):
                (row, col) = rowColToRequest[direction]
                nbrBinnedData = BinnedData(nbrFilename, row=row, col=col)
                binnedData.mergeNeighbour(nbrBinnedData, direction)

    return binnedData


###################################################################################################
def laspy2rec(infile):
    """
    Convert a LAS/LAZ file to a numpy record array using laspy.

    Parameters:
        infile (str or Path): Path to the input LAS/LAZ file.

    Returns:
        tuple: A tuple containing:
            - las_data (np.recarray): A structured numpy record array with the following fields:
                - "RETURN_NUMBER": Return number of the pulse (u1).
                - "NUMBER_OF_RETURNS": Total number of returns for the pulse (u1).
                - "TIMESTAMP": GPS timestamp of the pulse (<f8).
                - "INTENSITY": Intensity of the return (<i4).
                - "CLASSIFICATION": Classification of the point (u1).
                - "X": X coordinate (<f8).
                - "Y": Y coordinate (<f8).
                - "Z": Z coordinate (<f8).
            - header (laspy.LasHeader): The LAS/LAZ file header.

    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    # Check if the input file exists
    if not Path(infile).is_file():
        raise FileNotFoundError(f"The file '{infile}' does not exist.")

    # Open the LAS/LAZ file using laspy
    las = laspy.read(infile)

    # Extract the header
    header = las.header

    # Convert LAS data to a structured numpy record array
    las_data = np.rec.fromarrays(
        [
            las.return_num,  # Return number
            las.num_returns,  # Number of returns
            las.gps_time,  # GPS timestamp
            las.intensity,  # Intensity
            las.classification,  # Classification
            las.x,  # X coordinate
            las.y,  # Y coordinate
            las.z,  # Z coordinate
        ],
        names=[
            "RETURN_NUMBER",
            "NUMBER_OF_RETURNS",
            "TIMESTAMP",
            "INTENSITY",
            "CLASSIFICATION",
            "X",
            "Y",
            "Z",
        ],
        formats=["u1", "u1", "<f8", "<i4", "u1", "<f8", "<f8", "<f8"],
    )

    # Sort the record array by the "TIMESTAMP" field
    las_data.sort(order="TIMESTAMP")

    # Clean up
    del las

    return las_data, header


###################################################################################################
def standardise_lasf(what, when, utmZone, stageCode, projectName, outdir,
        data, easting, northing, tile_s, out_tile_s, binSize, filename_Parent,
        classesToExclude, minZ, maxZ, skipexisting):
    """
    Using laspy, rename file using naming convention, add index, remove noise and write out
    supplied files to .laz

    Parameters:
        what (str): The 'what' field of output filename
        when (str): The 'when' field of output filename
        utmZone (int): UTM zone number of projection for output files
        stageCode (str): 3-char stage code for output files
        projectName (str): 6-char name of data project
        outdir (str): Output directory for processed files.
        data (laspy.LasData): LAS/LAZ data to process.
        easting (float): Easting coordinate of top-left corner of the tile.
        northing (float): Northing coordinate of top-left corner of the tile.
        tile_s (int): Size of the input tile (metres).
        out_tile_s (int): Size of the output tile (metres).
        binSize (float): Bin size for indexing (metres).
        filename_Parent (str): Parent filename for metadata tracking.
        classesToExclude (list): List of integer point classification values to exclude
                                 from the data
        minZ, maxZ (float): Min and max acceptable values (metres) for point height.
                            Heights outside this range will be discarded as errors.
        skipexisting (bool): If true, skip file if output file already exists

    """
    segments = np.arange(0, int(tile_s), int(out_tile_s))

    projectionCode = qvf.makeProjectionCode(utmZone)
    stageAndZone = f"{stageCode}{projectionCode}"
    outfileTemplate = qvf.assemblefields([what, 'TILENAME', when, stageAndZone])
    outfileTemplate = qvf.setoptionfield(outfileTemplate, 'p', projectName)
    outfileTemplate = qvf.setsuffix(outfileTemplate, 'laz')
    outfileTemplate = os.path.join(outdir, outfileTemplate)

    for tile_x in segments:
        for tile_y in segments:
            northing_new = int(northing - tile_y)
            easting_new = int(easting + tile_x)
            tileWhere = qvf.makeTileWhere(easting, northing, utmZone)
            outFile = qvf.setwhere(outfileTemplate, tileWhere)

            if not (skipexisting and os.path.exists(outFile)):
                # Exclude points which are outside the tile to be output
                good_indices = (
                    (
                        (float(easting_new + out_tile_s - 0.001) > data.x)
                        & (float(easting_new + 0.001) <= data.x)
                    )
                    & (
                        ((northing_new - 0.001) >= data.y)
                        & ((northing_new - out_tile_s + 0.001) < data.y)
                    )
                )
                # Exclude points outside acceptable height range
                good_indices = (good_indices & ((minZ < data.z) & (data.z < maxZ)))
                # Exclude any point classes the user requested
                for classVal in classesToExclude:
                    good_indices = (good_indices & (data.classification != classVal))

                if np.sum(good_indices) > 0:
                    data2 = data[good_indices]
                    new_hdr = copy(data.header)
                    new_hdr.point_count = 0
                    new_las = laspy.LasData(new_hdr)

                    ##########################################################################
                    # GENERATE INDEX
                    nbinsRow = np.round(out_tile_s / binSize)  # ****
                    xIdx = ((data2.x - easting_new) // binSize).astype(np.int32)
                    yIdx = ((northing_new - np.array(data2.y)) // binSize).astype(
                        np.int32
                    )
                    index = ((yIdx * nbinsRow) + xIdx).astype(int)
                    nbins = int(nbinsRow * nbinsRow)

                    # Sort the points to group all points by bin. Create an index recording
                    # start and end points for each bin, allowing fast retrieval of points by bin.
                    # Note that we preserve the order of points within each bin.
                    sorter = np.argsort(index, stable=True)
                    data2 = data2[sorter]
                    binCounts = np.bincount(index, minlength=nbins)
                    bounds = np.cumsum(binCounts)
                    newIdx = np.concatenate(([0], bounds))

                    ##########################################################################
                    newIdx = np.array(newIdx, dtype=np.uint64)
                    nElems = int(nbins + 1)
                    #############################
                    binSizePacked = struct.pack("<dQ", np.float64(binSize), np.uint64(nbins))
                    # Instantiate and store new VLR.
                    new_vlr = laspy.VLR(
                        user_id=const.VLR_USERID_JRSRP,
                        record_id=const.VLR_RECORDID_BINSIZE,
                        description="Bin size & count",
                        record_data=binSizePacked,
                    )
                    new_las.vlrs.append(new_vlr)
                    #############################
                    binBoundsPacked = struct.pack(f"<{nElems}Q", *newIdx[0:nElems])
                    # Instantiate and store new VLR.
                    new_vlr = laspy.VLR(
                        user_id=const.VLR_USERID_JRSRP,
                        record_id=const.VLR_RECORDID_BINBOUNDS,
                        description="Bin start/end values",
                        record_data=binBoundsPacked,
                    )
                    new_las.vlrs.append(new_vlr)
                    #####################################
                    # check nbins match tile_s
                    test_bins = int((out_tile_s / binSize) ** 2)

                    if test_bins != nbins:
                        msg = "mis-match in bin indexing"
                        logger.error(msg)
                        raise ValueError(msg)
                    ###################################
                    new_las.points = data2.points.copy()
                    new_las.write(outFile)

                    # check file is not corrupted.
                    with laspy.open(outFile) as las:
                        point_count = las.header.point_count
                        xmin, ymin, zmin = las.header.min  # noqa
                        xmax, ymax, zmax = las.header.max  # noqa

                    if point_count < 1:
                        msg = f"File likely experienced wrapping: {outFile}"
                        logger.error(msg)
                        raise ValueError(msg)
                    if xmin < 1:
                        msg = f"File likely experienced wrapping: {outFile}"
                        logger.error(msg)
                        raise ValueError(msg)
                    if ymax < 1:
                        msg = f"File likely experienced wrapping: {outFile}"
                        logger.error(msg)
                        raise ValueError(msg)
                else:
                    data2 = None
            else:
                print(f"Skipping outfile {outFile}, as it already exists")

    del data


###################################################################################################
# LAZ INDEX
def read_laz_index(infile, tile_s):
    """
    e.g binSize,nbins,newIdx,las_data = read_laz_index(outfn)

    Reads the binning index from a LAS/LAZ file.

    Returns:
    - binSize: The size of each bin.
    - nbins: The number of bins.
    - newIdx: An array of indices indicating the bin positions.
    - las_data: The LAS data as a structured NumPy array.

    """

    if (Path(infile)).is_file():
        las = laspy.read(infile)
        las_data = np.rec.fromarrays(
            [
                las.return_num,
                las.num_returns,
                las.gps_time,
                las.intensity,
                las.classification,
                las.x,
                las.y,
                las.z,
            ],
            names=[
                "RETURN_NUMBER",
                "NUMBER_OF_RETURNS",
                "TIMESTAMP",
                "INTENSITY",
                "CLASSIFICATION",
                "X",
                "Y",
                "Z",
            ],
            formats=["u1", "u1", "<f8", "<i4", "u1", "<f8", "<f8", "<f8"],
        )

        inVLRs = las.vlrs
        # Initialize variables
        binSize = None
        nbins = None
        newIdx = None

        # Search for the required VLRs
        bin_info_vlr = None
        bin_index_vlr = None

        #########################################################################
        # Look for the conformant VLRs
        bin_info_vlr_list = inVLRs.get_by_id(user_id=const.VLR_USERID_JRSRP,
                                             record_ids=[const.VLR_RECORDID_BINSIZE])
        if len(bin_info_vlr_list) > 0:
            # These are the newer, standard-conformant VLRs for the index
            bin_info_vlr = bin_info_vlr_list[0]
            bin_index_vlr_list = inVLRs.get_by_id(user_id=const.VLR_USERID_JRSRP,
                                                  record_ids=[const.VLR_RECORDID_BINBOUNDS])
            bin_index_vlr = bin_index_vlr_list[0]
            oldVlrs = False
        else:
            # Old non-conformant VLRs
            for vlr in inVLRs:
                if vlr.user_id == "BINSIZE":
                    # VLR containing binSize and nbins
                    bin_info_vlr = vlr
                elif vlr.user_id == "BIN_POS":
                    # VLR containing bin positions
                    bin_index_vlr = vlr
            oldVlrs = True
        #########################################################################

        # Check if required VLRs were found
        if bin_info_vlr is None:
            raise ValueError(
                f"Could not find VLR for bin size & count in file {infile}"
            )

        if bin_index_vlr is None:
            raise ValueError(
                f"Could not find VLR for bin bounds in file {infile}"
            )

        # Unpack binSize and nbins
        bin_info_format = "<dQ"
        if oldVlrs:
            bin_info_format = "2d"
        binSize, nbins = struct.unpack(bin_info_format, bin_info_vlr.record_data)

        if nbins > 1:
            n_elements = int(nbins + 1)
            bin_index_format = f"<{n_elements}Q"
            if oldVlrs:
                bin_index_format = f"{n_elements}d"
            newIdx = struct.unpack(bin_index_format, bin_index_vlr.record_data)
            newIdx = np.array(newIdx, dtype=np.int32)
        else:
            msg = f"Failed to import laz binning index in file {infile}"
            raise ValueError(msg)
        del las

        return binSize, nbins, newIdx, las_data
    else:
        msg = f"File {infile} not found"
        raise FileNotFoundError(msg)


###################################################################################################
def bin_data(infile, tile_s, nbins):
    """
    Process LAS/LAZ files by binning data into a grid structure.

    Parameters:
        infile (str): Name of the input LAS/LAZ file to process.
        tile_s (int): Tile size.
        nbins (int): Number of bins to divide the tile into.

    Returns:
        tuple: A tuple containing:
            - bins (dict): A dictionary where keys are bin names (e.g., "row_<row>_col_<col>")
                           and values are data arrays.
            - yst (int): Starting Y-coordinate of the tile.
            - xst (int): Starting X-coordinate of the tile.
            - rowS (int): Number of rows in the grid.
    """
    # Get bin indices and related metadata
    tile_idx, bins2access, xbinID, ybinID, rowS = get_bin_indices(nbins, tile_s)

    # Create bin names for the grid
    bin_names = [
        f"row_{row}_col_{col}" for row in range(rowS + 2) for col in range(rowS + 2)
    ]
    bin_names = np.array(bin_names)
    bins = {pos: np.zeros(1) for pos in bin_names}

    # Parse metadata from the input file name
    where = qvf.getwhere(infile)
    (easting, northing, utmZone) = filenaming_methods.decomposeWhereField(where)

    # Process neighbouring tiles
    for pos, p in enumerate(range(8)):
        # Determine the location of the neighbouring tile
        (xoffset, yoffset) = tile_idx[p]
        neighbourWhere = filenaming_methods.neighbourTileWhere(where, xoffset, yoffset)
        neighbourfile = qvf.setwhere(infile, neighbourWhere)
        # If the neighbouring file exists, process it
        if Path(neighbourfile).is_file():
            binSize, nbins, newIdx, data = read_laz_index(neighbourfile, tile_s)
            bb = bins2access[pos]

            # Assign data to the appropriate bins
            for loc, val in enumerate(bb):
                bin_name = f"row_{ybinID[pos][loc]}_col_{xbinID[pos][loc]}"
                bins[bin_name] = np.copy(data[newIdx[int(val)] : newIdx[int(val + 1)]])
            del data

    # Process the main input file
    binSize, nbins, newIdx, data = read_laz_index(infile, tile_s)  # noqa

    # Assign data to bins for the main tile
    ct = 0
    for row in np.flip(range(1, rowS + 1)):
        for col in range(1, rowS + 1):
            bin_name = f"row_{row}_col_{col}"
            bins[bin_name] = np.copy(data[newIdx[ct] : newIdx[ct + 1]])
            ct += 1
    del data

    return bins, northing, easting, rowS


###################################################################################################
def get_bin_indices(nbins, tile_s):
    """
    Generate indices for bins of surrounding eight LAS/LAZ tiles for processing.

    Parameters:
        nbins (int): Total number of bins in the tile.
        tile_s (int): Tile size.

    Returns:
        tuple: A tuple containing:
            - tile_idx (np.ndarray): Array of offsets for the eight surrounding tiles.
            - bins2access (list): List of bin indices to access for each surrounding tile.
            - xbinID (list): List of X bin IDs for each surrounding tile.
            - ybinID (list): List of Y bin IDs for each surrounding tile.
            - rowS (int): Number of rows (or columns) in the grid.
    """
    # Define offsets for the eight surrounding tiles
    tile_idx = (
        np.array([[-1, 1], [0, 1], [1, 1], [-1, 0], [1, 0], [-1, -1], [0, -1], [1, -1]])
        * tile_s
    )

    # Calculate the number of rows (or columns) in the grid
    rowS = int(np.sqrt(nbins))

    # Define bins to be accessed for each surrounding tile
    bins2access = [
        [nbins - 1],  # Top-right corner
        np.arange((nbins - rowS), nbins),  # Top row
        [(nbins - rowS)],  # Top-left corner
        np.arange((rowS - 1), rowS**2, rowS),  # Right column
        np.arange(rowS) * rowS,  # Left column
        [rowS - 1],  # Bottom-right corner
        np.arange(rowS),  # Bottom row
        [0],  # Bottom-left corner
    ]

    # Define X bin IDs for each surrounding tile
    xbinID = [
        [0],  # Top-right corner
        np.arange(1, rowS + 1),  # Top row
        [rowS + 1],  # Top-left corner
        np.zeros(rowS, dtype=np.int32),  # Right column
        np.zeros(rowS, dtype=np.int32) + rowS + 1,  # Left column
        [0],  # Bottom-right corner
        np.arange(1, rowS + 1),  # Bottom row
        [rowS + 1],  # Bottom-left corner
    ]

    # Define Y bin IDs for each surrounding tile
    ybinID = [
        [rowS + 1],  # Top-right corner
        np.zeros(rowS, dtype=np.int32) + rowS + 1,  # Top row
        [rowS + 1],  # Top-left corner
        np.flip(np.arange(1, rowS + 1)),  # Right column
        np.flip(np.arange(1, rowS + 1)),  # Left column
        [0],  # Bottom-right corner
        np.zeros(rowS, dtype=np.int32),  # Bottom row
        [0],  # Bottom-left corner
    ]

    return tile_idx, bins2access, xbinID, ybinID, rowS


###################################################################################################
# def add_hag_evlr(hag_data,outf):
#     """
#     code for adding height-above-ground (hag) to evlr

#     can store as int with only cm precision required

#     just dumped example for now
#     """
#     a = [i for i in hag_data][:]
#     hag_pk = struct.pack(f"{len(hag_data)}i",*a)

#     new_evlr = laspy.vlrs.vlr.VLR(
#                 user_id="HeightAboveGround",   # Max 16 characters
#                 record_id=12345,            # Unsigned short integer ID
#                 description="Custom Metadata", # Max 32 characters
#                 record_data=hag_pk # Must be raw bytes
#                 )

#     if outlas.evlrs is None:
#         outlas.evlrs = laspy.vlrs.vlrlist.VLRList()

#     outlas.evlrs.append(new_evlr)
#     outlas.write(outf)

# def read_evlr(fn):
#     """
#     test reading evlr
#     just dumped example for now
#     """
#     las = laspy.read(fn)
#     for evlr in las.evlrs:
#         print(evlr.user_id, evlr.record_id)
#         byte_data = list(evlr.record_data)
#         print(len(byte_data))
#         print(byte_data[0:9])


# ##############################################################################################################################

# def pdal_convert_epsg():
#     """
#     just a placeholder for now
#     """
#     a=a

# ##############################################################################################################################

# def check_las_version(fn,outf,point_format_id=6):
#     """
#     check las verion - if not version 1.4 then convert
#     """
#     las_data = laspy.read(fn)
#     new_las = laspy.convert(las_data, point_format_id=point_format_id, file_version="1.4")
#     new_las.write(outf)


#########################################
##
def run_zipf(fn, laz_flist, outdr, stagecode="ba2"):
    """
    moved here - needs to be updated
    Create zip files for BA2 and BA3 products.

    Parameters:
        fn: laz file named following name convention
        laz_flist: full path to list of laz files
        outdr: directory to write zip file to
        stagecode:

    """
    # Read the list of LAS/LAZ files
    with open(laz_flist) as f:
        fns = [line.strip() for line in f]

    # Change the working directory to the directory containing the laz_flist
    dr = Path(laz_flist).parent
    os.chdir(dr)

    # Validate the provided filename
    if not fn:
        print("Error: --fn argument is required.")
        sys.exit(1)

    # tileBasename = Path(fn).with_suffix('')
    tileDict = filenaming_methods.createTileDict(fn, 1000)

    # Generate filenames for BA2 and BA3 zip files
    what = f"ap{tileDict['instrument']}{tileDict['returntype']}"
    where = f"r{tileDict['project']}"
    when = f"{tileDict['date']}"
    stageAndZone = f"{stagecode}{tileDict['zone_prefix']}{tileDict['zoneCode']}"
    fn_out = f"{what}_{where}_{when}_{stageAndZone}.zip"

    # Write zip file
    archive_name = Path(outdr).joinpath(fn_out)
    with zipfile.ZipFile(archive_name, "w") as zf:
        for lazfile in fns:
            zf.write(lazfile, arcname=Path(lazfile).name)
    print(f"{fn_out} zip file created: {archive_name}")
