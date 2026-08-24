#!/usr/bin/env python

"""
management of file naming - but probably will be superceded?

File Naming Structure:
    <what>_<where>_<when>_<processing>[_optional][.suffix]

What:
    - `ssiipp` where:
        - `ss`: Platform type
            - 'ap': Airborne platform
            - 'gp': Ground platform
            - 'sb': Spaceborne
            - 'is': IceSAT
        - `ii`: Sensor code (see `get_sensorCodes` for details)
        - `pp`: Product type
            - `dr`: Discrete return
            - `wf`: Waveform
            - `cw`: Continuous wave (phase-shift) LiDAR

Where:
    - Tile format: `x<easting>ys<northing>z<zone_number>` (e.g., `x561000ys7155000z55`)
    - Mosaic filename: Six-character region name prefixed with 'r'.

When:
    - For ALS captures, the time refers to the year(s) of capture.

Examples:
    Tile:
        `apr8dr_x448750ys7133000z56_2021_ba1m6_pmivasf.laz`
        - `ap`: Airborne platform
        - `r8`: Sensor code
        - `dr`: Discrete return

    Output Product:
        `apr8dr_rmivasf_2021_bbpm6_r50cm.tif`
"""

# Configure logging
import logging
import os

import numpy as np
from osgeo import gdal

from airborne_ls import qvf

logger = logging.getLogger(__name__)


stageByProductName = {
    "dem" : "bb0",
    "maxH" : "bb1",
    "intens" : "bb2",
    "grdR" : "bb3",
    "NonGrdCodes" : "bb4",     # Non-ground codes
    "fstDens" : "bb5",         # First return density
    "percentile1" : "bb8",
    "percentile5" : "bb9",
    "percentile25" : "bba",
    "percentile50" : "bbb",
    "percentile75" : "bbc",
    "percentile95" : "bbd",
    "percentile99" : "bbe",
    "fpc" : "bbh",             # Foliage profile curve
    "demHS" : "bbi",           # DEM hillshade
    "csm" : "bbm",             # Canopy surface model
    "chm" : "bbn",             # Canopy height model
}
# Reserved for possible future inclusion
#    "demFilled4hydro": "bbp"   # DEM filled for hydrological analysis
#    "flowAccumulation": "bbq"
#    "DFMEraw": "bbr"           # Digital Flow Model Elevation (raw)
#    "rgDFME": "bbs"            # Region-growing DFME
#    "depressionDepth": "bbt"
#    "slopedh5x5": "bbu"        # Slope derived from 5x5 window
#    "streamClassi": "bbv"      # Stream classification

# And a reverse lookup of the same information
productNameByStage = {stageByProductName[k]: k for k in stageByProductName}


def get_stageDict():
    """
    Returns a dictionary mapping stage codes to product descriptions.

    Returns:
        dict: A dictionary where keys are stage codes (e.g., 'bb0') and
              values are product descriptions.
    """
    productDict = {
        "bb0": "dem",
        "bb1": "maxH",
        "bb2": "intens",
        "bb3": "grdR",
        "bb4": "NonGrdCodes",  # Non-ground codes
        "bb5": "fstDens",  # First return density
        "bb8": "percentile1",
        "bb9": "percentile5",
        "bba": "percentile25",
        "bbb": "percentile50",
        "bbc": "percentile75",
        "bbd": "percentile95",
        "bbe": "percentile99",
        "bbh": "fpc",  # Foliage profile curve
        "bbi": "demHS",  # DEM hillshade
        "bbm": "csm",  # Canopy surface model
        "bbn": "chm",  # Canopy height model
        # "bbp": "dem_filled4hydro",  # DEM filled for hydrological analysis
        # "bbq": "flowAccumulation",
        # "bbr": "DFME_raw",  # Digital Flow Model Elevation (raw)
        # "bbs": "rg_DFME",  # Region-growing DFME
        # "bbt": "depression_Depth",
        # "bbu": "slope_dh_5x5",  # Slope derived from 5x5 window
        # "bbv": "streamClassi",  # Stream classification
    }

    return productDict


def get_recordID_Dict():
    recordID_Dict = {
        "run_lidar_standardisation": 998,
        # "???": 999,
        "bb0": 1000,
        "bb1": 1001,
        "bb2": 1002,
        "bb3": 1003,
        "bb4": 1004,
        "bb5": 1005,
        "bb8": 1006,
        "bb9": 1007,
        "bba": 1008,
        "bbb": 1009,
        "bbc": 1010,
        "bbd": 1011,
        "bbe": 1012,
        "bbh": 1013,
        "bbi": 1014,
        "bbm": 1015,
        # "bbp": 1016,
        # "bbq": 1017,
        # "bbr": 1018,
        # "bbs": 1019,
        # "bbt": 1020,
        # "bbu": 1021,
        # "bbv": 1022,
        "bbn": 2024,
    }

    return recordID_Dict


def get_psizeDict(psize=0.5, ptile_s=5, fpc_psize=10, chm_psize=0.2):
    """
    Generate a dictionary mapping product codes to their respective resolutions.

    Parameters:
        psize (float): Resolution for DEM and maxH in metres (default: 0.5).
        ptile_s (float): Resolution for percentiles in metres (default: 5).
        fpc_psize (float): Resolution for FPC in metres (default: 10).
        chm_psize (float): Resolution for CHM in metres (default: 0.2).

    Returns:
        dict: A dictionary where keys are product codes (e.g., 'bb0') and
              values are resolution strings.
    """
    # Convert resolutions to appropriate string formats
    psize = f"r{int(psize * 100)}cm"
    chm_psize = f"r{int(chm_psize * 100)}cm"
    if ptile_s < 10:
        ptile_s = f"r{int(ptile_s * 100)}cm"
    else:
        ptile_s = f"r{int(ptile_s)}m"
    fpc_psize = f"r{int(fpc_psize)}m"

    # Print CHM resolution for debugging purposes
    print(chm_psize)

    # Create the dictionary mapping product codes to resolutions
    psizeDict = {
        "bb0": psize,  # DEM
        "bb1": psize,  # maxH
        "bb2": psize,  # intensity
        "bb3": psize,  # grdR
        "bb4": psize,  # Non-ground codes
        "bb5": psize,  # First return density
        "bb8": ptile_s,  # 1st percentile
        "bb9": ptile_s,  # 5th percentile
        "bba": ptile_s,  # 25th percentile
        "bbb": ptile_s,  # 50th percentile
        "bbc": ptile_s,  # 75th percentile
        "bbd": ptile_s,  # 95th percentile
        "bbe": ptile_s,  # 99th percentile
        "bbg": psize,  # Additional DEM product
        "bbh": fpc_psize,  # Foliage profile curve (FPC)
        "bbi": psize,  # DEM hillshade
        "bbj": psize,  # Additional product
        "bbm": psize,  # Canopy surface model (CSM)
        "bbn": chm_psize,  # Canopy height model (CHM)
        # "bbo": psize,  # Additional product
        # "bbp": psize,  # DEM filled for hydrological analysis
        # "bbq": psize,  # Flow accumulation
        # "bbr": psize,  # Digital Flow Model Elevation (raw)
        # "bbs": psize,  # Region-growing DFME
        # "bbt": psize,  # Depression depth
        # "bbu": psize,  # Slope derived from 5x5 window
        # "bb6": psize,  # Additional product
        # "bb7": psize,  # Additional product
        # "bbf": psize,  # Additional product
        # "bbv": psize,  # Stream classification
    }

    return psizeDict


def resolutionStrFromMetres(metres):
    """
    Create a resolution string from the given number of metres
    """
    if metres < 1:
        cm = int(metres * 100)
        resStr = f"{cm}cm"
    else:
        m = int(metres)
        resStr = f"{m}m"
    return resStr


def decomposeWhereField(where):
    """
    Break up the given where field and return the values of the sub-fields

    Parameters:
      where: A where field, of the form x<xCoord>ys<yCoord>z<utmZone>
             Note that the 's' after the 'y' indicate UTM South

    Returns:
      (x, y, zone): The (X, Y) coordinates (metres) and the UTM zone (negative for South)
    """
    xNdx = where.find('x')
    yNdx = where.find('y')
    zNdx = where.find('z')
    xCoord = int(where[xNdx + 1:yNdx])
    yCoord = int(where[yNdx + 2:zNdx])
    utmZone = int(where[zNdx + 1:])
    if where[yNdx + 1] == 's':
        utmZone = -utmZone
    return (xCoord, yCoord, utmZone)


def neighbourTileWhere(where, xOffset, yOffset):
    """
    Given the where field of a file name for a single tile of data, return the
    where field for a neighbouring tile, based on the xOffset & yOffset parameters.

    The X & Y offsets are taken to be in metres, and are exactly one tile in
    some direction. So, for example, if the tile size is 1000m, then xOffset of -1000
    would indicate the tile to the west, and a yOffset of +1000 would indicate the
    tile to the north.

    Parameters:
      where (str): Where field for a single tile
      xOffset (int): Offset (metres) in X direction to top-left of neighbour tile
      yOffset (int): Offset (metres) in Y direction to top-left of neighbour tile

    Returns:
      nbrwhere (str): Where field of requested neighbouring tile
    """
    (xCoord, yCoord, utmZone) = decomposeWhereField(where)
    newX = xCoord + int(xOffset)
    newY = yCoord + int(yOffset)
    newWhere = qvf.makeTileWhere(newX, newY, utmZone)
    return newWhere


def get_outfnames(outputBasename, psize=0.5, ptile_s=5, fpc_psize=10, chm_psize=0.2,
        pptiles=(1, 5, 25, 50, 75, 95, 99), driverName='GTiff'):
    """
    Generate filenames for individually processed tiles and intermediate products.

    Parameters:
        outputBasename (str): Base name for the output files.
        psize (float): Resolution for DEM, maxH, and other products in metres (default: 0.5).
        ptile_s (float): Resolution for percentile tiles in metres (default: 5).
        fpc_psize (float): Resolution for FPC in metres (default: 10).
        chm_psize (float): Resolution for CHM in metres (default: 0.2).
        pptiles (list[int]): List of percentiles to generate filenames for
                             (default: (1, 5, 25, 50, 75, 95, 99)).

    Returns:
        dict: A dictionary where keys are product names and
              values are their corresponding filenames.
    """
    stdRes = resolutionStrFromMetres(psize)
    fpcRes = resolutionStrFromMetres(fpc_psize)
    chmRes = resolutionStrFromMetres(chm_psize)
    pcntileRes = resolutionStrFromMetres(ptile_s)

    fnames = {}

    for productName in stageByProductName:
        stageCode = stageByProductName[productName]
        outfile = qvf.setstagecode(outputBasename, stageCode)
        outfile = qvf.setoptionfield(outfile, 'l', productName)

        # Work out which resolution string to use
        resStr = stdRes
        if productName == "fpc":
            resStr = fpcRes
        elif productName == "chm":
            resStr = chmRes
        elif productName.startswith('percentile'):
            resStr = pcntileRes

        outfile = qvf.setoptionfield(outfile, 'r', resStr)
        suffix = getSuffixFromDriverName(driverName)
        outfile = qvf.setsuffix(outfile, suffix)
        fnames[productName] = outfile

    # Return the filenames dictionary
    return fnames


def getSuffixFromDriverName(driverName):
    """
    Get the preferred suffix for the given GDAL driver name
    """
    drvr = gdal.GetDriverByName(driverName)
    suffix = drvr.GetMetadataItem('DMD_EXTENSION')
    return suffix


def createTileDict(tile, tile_s):
    """
    Parse the name of an input tile to extract metadata for output products.

    Example:
        Input tile: apr8dr_x448750ys7133000_2021_ba1m6_pmivasf.laz

    Parameters:
        tile (str): The filename of the input tile.
        tile_s (int): The size of the tile.

    Returns:
        dict: A dictionary containing parsed metadata from the tile name.
    """
    # Split the tile name into components
    components = tile.split("_")
    if len(components) != 5:
        msg = f"Input {tile} error in filename components."
        logger.error(msg)
        raise ValueError(msg)

    # Extract components from the filename
    what = components[0]
    where = components[1]
    when = components[2]
    processing = components[3]

    # Validate the 'what' component
    if what[:2] != "ap":
        msg = f"Filename sensor check failed.Input {tile} error in filename components."
        logger.error(msg)
        raise ValueError(msg)

    # Parse the project name
    project = components[4].split(".")[0]
    if len(project[1:]) != 6:
        msg = (
            f"Invalid project name '{project}'. "
            f"Expected 1 prefix character followed by 6 characters."
        )
        logger.error(msg)
        raise ValueError(msg)

    # Create the tile dictionary
    tileDict = {
        "rawtilename": os.path.basename(tile),  # Original tile filename
        "satellite": what[:2],  # Satellite/platform code (e.g., 'ap')
        "instrument": what[2:4],  # Instrument code
        "returntype": what[4:6],  # Return type (e.g., 'dr')
        "date": int(when),  # Year of capture
        "tile": where,  # Tile location (e.g., x448750ys7133000)
        "project": project[1:],  # Project code (6 characters)
        "stage": processing[:3],  # Processing stage code
        "zone_prefix": processing[-2:-1],  # Zone prefix
        "zoneCode": int(processing[-1:]),  # Zone code
        "xst": np.int32(where[1:7]),  # Starting x-coordinate
        "yst": np.int32(where[9:16]),  # Starting y-coordinate
        "tile_s": np.int32(tile_s),  # Tile size
        "components": components,  # Full list of filename components
    }

    return tileDict
