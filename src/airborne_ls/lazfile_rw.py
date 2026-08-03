#!/usr/bin/env python

"""
las/laz file manipulations for lidarveg processing of airborne lidar

key activites / options: 

-- check file for issues with supplied file
-- build a processing index
-- standardise laz version converting from las where required 
-- update projections where needed 

"""

import struct
import sys
from pathlib import Path
from copy import copy

import laspy
import numpy as np

################################################################################################################
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
            las.return_num,        # Return number
            las.num_returns,       # Number of returns
            las.gps_time,          # GPS timestamp
            las.intensity,         # Intensity
            las.classification,    # Classification
            las.x,                 # X coordinate
            las.y,                 # Y coordinate
            las.z,                 # Z coordinate
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
    

########################################################################################################################################
def standardise_lasf(fn_base,outdr,data, easting, northing, tile_s, out_tile_s, binSize, filename_Parent):
    """
        Using laspy, rename file using naming convention, add index, remove noise and write out supplied files to .laz

        Parameters:
            fn_base (str): Base filename for output files.
            outdr (str): Output directory for processed files.
            data (laspy.LasData): LAS/LAZ data to process.
            easting (float): Easting coordinate of the tile.
            northing (float): Northing coordinate of the tile.
            tile_s (float): Size of the input tile (metres).
            out_tile_s (float): Size of the output tile (metres).
            binSize (float): Bin size for indexing (metres).
            filename_Parent (str): Parent filename for metadata tracking.
            
        Returns:
            str: Status message indicating the result of the processing.
    """
    header = data.header       
    test_divisible = tile_s % out_tile_s == 0
    
    if test_divisible:
        segments=(np.arange(0,int(tile_s),int(out_tile_s)))
        
        for tile_x in segments:
            for tile_y in segments:
    
                northing_new = int(northing-tile_y)
                easting_new = int(easting+tile_x)
        
                pts = fn_base.split('_')
                fn_where = "x%sys%s" % (easting_new, northing_new)
                fn_base2 = ("_").join([pts[0], fn_where, pts[2], pts[3], pts[4]])                
                outfn = str(Path(outdr).joinpath(fn_base2))
                                
                # simply exclude points outside the tile extents + irrelevant codes/data
                # 7 = low point noise, 18 = high point noise.. note some providers can use different / new codes
                ## z-thresholds are problematic.. as you can have negative elevation and what upper limit? for aus 3000m works 

                good_indices =  ( (float(easting_new+out_tile_s-.001) > data.x)  & (float(easting_new+.001) <=  data.x) ) &\
                                ( ((northing_new-.001) >= data.y) & ((northing_new-out_tile_s+.001) < data.y) )  &\
                                ( (float(-10.) < data.z) &  (float(3000.) > data.z) )  &\
                                  (data.classification != 7)  & (data.classification != 18) & (data.classification !=64) 
                
                if np.sum(good_indices)>0:                                                            
                    
                    data2 = data[good_indices]                                    
                    new_hdr = copy(data.header)
                    new_hdr.point_count = 0
                    new_las = laspy.LasData(new_hdr)
                    
                    ##########################################################################
                    ## GENERATE INDEX
                    nbinsRow = np.round(out_tile_s / binSize) ## ****
                    xIdx = ((data2.x - easting_new) // binSize).astype(np.int32)
                    yIdx = ((northing_new - np.array(data2.y)) // binSize).astype(np.int32)
                    index = ((yIdx * nbinsRow) + xIdx).astype(int)
                    nbins = int(nbinsRow * nbinsRow)
                    
                    start = 0
                    newLaz = []
                    newIdx = [0]
                    sortingIdx = []
                    ### needs speeding up..
                    for ct, binIdx in enumerate(range(nbins + 1)):
                        vals2 = np.argwhere(index == binIdx)
                        if sum(vals2) >= 0:
                            start += len(vals2)
                            if len(sortingIdx) < 1:
                                sortingIdx = vals2
                            else:
                                sortingIdx = np.append(sortingIdx, vals2)
                            newIdx.append(np.copy(start))
                    newIdx.append(len(data2.x))    
                    data2 = data2[sortingIdx]   
                    del sortingIdx 
                    
                    ##########################################################################
                    newIdx = np.array(newIdx, dtype=np.float64)                      
                    nElems = int(nbins + 1)
                    #############################
                    binS = struct.pack("2d", np.float64(binSize), np.float64(nbins))
                    # Instantiate and store new VLR.
                    new_vlr = laspy.VLR(user_id="BINSIZE", record_id=1, \
                                        description="defines binSize chunk", record_data=binS)                    
                    new_las.vlrs.append(new_vlr)
                    #############################
                    bin_pos = struct.pack("%sd" % (nElems), *newIdx[0:nElems])
                    # Instantiate and store new VLR.
                    new_vlr = laspy.VLR(user_id="BIN_POS", record_id=2, \
                                        description="defines BIN_POS", record_data=bin_pos)                    
                    new_las.vlrs.append(new_vlr)
                    #####################################
                    # check nbins match tile_s
                    test_bins = int((out_tile_s / binSize)**2)
                    try:
                        assert test_bins==nbins, "mis-match in bin indexing"
                    except AssertionError as err:
                        raise err
                    ###################################
                    new_las.points = data2.points.copy()
                    new_las.write(outfn)          
       
                    # check file is not corrupted.                   
                    with laspy.open(outfn) as las:    
                        point_count = las.header.point_count
                        xmin,ymin,zmin = las.header.min
                        xmax,ymax,zmax = las.header.max
                        
                    try:
                        assert point_count > 0, "File likely experienced wrapping: %s" % (outfn)
                    except AssertionError as err:
                        raise err
                    try:
                        assert xmin > 1, "File likely experienced wrapping: %s" % (outfn)
                    except AssertionError as err:
                        raise err
                    try:
                        assert ymax > 1, "File likely experienced wrapping: %s" % (outfn)
                    except AssertionError as err:
                        raise err
                    status = 'file indexed'
                else:                    
                    data2 = None
                    status = "Status: No Good Indices"
    else:
        print('tile dimensions are not divisible %s and %s') % (tile_s, out_tile_s)
        
    del data
    
    return status

##############################################################################################################################
## LAZ INDEX 
def read_laz_index(infile,tile_s):
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
        header = las.header
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
        for vlr in inVLRs:
            if vlr.user_id == 'BINSIZE':
                # VLR containing binSize and nbins
                bin_info_vlr = vlr
            elif vlr.user_id == 'BIN_POS':
                # VLR containing bin positions
                bin_index_vlr = vlr
        #########################################################################
        
        # Check if required VLRs were found
        if bin_info_vlr is None:
            raise ValueError(f"Could not find VLR with user_id='BINSIZE' and record_id=1 in file {infile}")

        if bin_index_vlr is None:
            raise ValueError(f"Could not find VLR with user_id='BIN_POS' and record_id=2 in file {infile}")

        # Unpack binSize and nbins
        binSize, nbins = struct.unpack("2d", bin_info_vlr.record_data)

        if nbins > 1:
            n_elements = int(nbins + 1)
            newIdx = struct.unpack(f"{n_elements}d", bin_index_vlr.record_data)
            newIdx = np.array(newIdx, dtype=np.int32)
        else:
            msg = f"Failed to import laz binning index in file {infile}"
            raise ValueError(msg)
        del las

        return binSize, nbins, newIdx, las_data
    else:        
        msg = "File %s not found" % (infile)
        raise FileNotFoundError(msg)

##############################################################################################################################
def bin_data(indir, infile, tile_s, nbins):
    """
    Process LAS/LAZ files by binning data into a grid structure.

    Parameters:
        indir (str): Directory containing the input LAS/LAZ files.
        infile (str): Name of the input LAS/LAZ file to process.
        tile_s (int): Tile size.
        nbins (int): Number of bins to divide the tile into.

    Returns:
        tuple: A tuple containing:
            - bins (dict): A dictionary where keys are bin names (e.g., "row_x_col_y") and values are data arrays.
            - yst (int): Starting Y-coordinate of the tile.
            - xst (int): Starting X-coordinate of the tile.
            - rowS (int): Number of rows in the grid.
    """
    # Get bin indices and related metadata
    tile_idx, bins2access, xbinID, ybinID, rowS = get_bin_indices(nbins, tile_s)

    # Create bin names for the grid
    bin_names = [
        f"row_{row}_col_{col}" for row in range(0, rowS + 2) for col in range(0, rowS + 2)
    ]
    bin_names = np.array(bin_names)
    bins = {pos: np.zeros((1)) for pos in bin_names}

    # Construct the full path to the input file
    infileFull = Path(indir).joinpath(infile)

    # Parse metadata from the input file name
    fn_dict = prod_stagecodes.createTileDict(infile, tile_s)

    # Process neighbouring tiles
    for pos, p in enumerate(range(8)):
        # Determine the location of the neighbouring tile
        where = f"x{int(fn_dict['xst'] + tile_idx[p, 0])}ys{int(fn_dict['yst'] + tile_idx[p, 1])}"
        fn = "_".join(
            (
                fn_dict["components"][0],
                where,
                fn_dict["components"][2],
                fn_dict["components"][3],
                fn_dict["components"][4],
            )
        )
        neighbourfile = os.path.join(indir, fn)

        # If the neighbouring file exists, process it
        if os.path.exists(neighbourfile):
            binSize, nbins, newIdx, data = read_laz_index(neighbourfile, tile_s)
            bb = bins2access[pos]

            # Assign data to the appropriate bins
            for loc, val in enumerate(bb):
                bin_name = f"row_{xbinID[pos][loc]}_col_{ybinID[pos][loc]}"
                bins[bin_name] = np.copy(data[newIdx[int(val)] : newIdx[int(val + 1)]])
            del data

    # Process the main input file
    binSize, nbins, newIdx, data = read_laz_index(infileFull, tile_s)

    # Assign data to bins for the main tile
    ct = 0
    for col in np.flip(range(1, rowS + 1)):
        for row in range(1, rowS + 1):
            bin_name = f"row_{row}_col_{col}"
            bins[bin_name] = np.copy(data[newIdx[ct] : newIdx[ct + 1]])
            ct += 1
    del data

    return bins, fn_dict["yst"], fn_dict["xst"], rowS

##############################################################################################################################
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

##############################################################################################################################
##############################################################################################################################
def add_hag_evlr(hag_data,outf):
    """
    code for adding height-above-ground (hag) to evlr  

    can store as int with only cm precision required 

    just dumped example for now
    """    
    a = [i for i in hag_data][:]
    hag_pk = struct.pack(f"{len(hag_data)}i",*a)

    new_evlr = laspy.vlrs.vlr.VLR(
                user_id="HeightAboveGround",   # Max 16 characters
                record_id=12345,            # Unsigned short integer ID
                description="Custom Metadata", # Max 32 characters
                record_data=hag_pk # Must be raw bytes
                )

    if outlas.evlrs is None:
        outlas.evlrs = laspy.vlrs.vlrlist.VLRList()

    outlas.evlrs.append(new_evlr)
    outlas.write(outf)

def read_evlr(fn):
    """
    test reading evlr
    just dumped example for now
    """
    las = laspy.read(fn)
    for evlr in las.evlrs:
        print(evlr.user_id, evlr.record_id)
        byte_data = list(evlr.record_data)
        print(len((byte_data)))
        print(byte_data[0:9])


##############################################################################################################################

def pdal_convert_epsg():
    """
    just a placeholder for now
    """
    a=a

##############################################################################################################################

def check_las_version(fn,outf,point_format_id=6):
    """
    check las verion - if not version 1.4 then convert
    """
    las_data = laspy.read(fn)
    new_las = laspy.convert(las_data, point_format_id=point_format_id, file_version="1.4")
    new_las.write(outf)



#########################################
##
def run_zipf(fn,laz_flist,outdr, stagecode='ba2'):
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
    fns = [line.strip() for line in open(laz_flist)]

    # Change the working directory to the directory containing the laz_flist
    dr = Path(laz_flist).parent
    os.chdir(dr)

    # Validate the provided filename
    if not fn:
        print("Error: --fn argument is required.")
        sys.exit(1)

    tileBasename = Path(fn).with_suffix('')
    tileDict = prod_stagecodes.createTileDict(fn, 1000)

    # Generate filenames for BA2 and BA3 zip files
    fn_out = f"ap{tileDict['instrument']}{tileDict['returntype']}_r{tileDict['project']}_{tileDict['date']}_{stagecode}{tileDict['zone_prefix']}{tileDict['zoneCode']}.zip"
    
    # Write zip file    
    archive_name = Path(outdr).joinpath(fn_out)
    with zipfile.ZipFile(archive_name, "w") as zf:
        for lazfile in fns:
            zf.write(lazfile, arcname=Path(lazfile).name)
    print(f"{fn_out} zip file created: {archive_name}")


    
