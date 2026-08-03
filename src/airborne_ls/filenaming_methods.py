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

logger = logging.getLogger(__name__)

def get_stageDict():
    """
    Returns a dictionary mapping stage codes to product descriptions.

    Returns:
        dict: A dictionary where keys are stage codes (e.g., 'bb0') and values are product descriptions.
    """
    productDict = {
        "bb0": "dem",
        "bb1": "maxH",
        "bb2": "intens",
        "bb3": "grdR",
        "bb4": "NonGrd_codes",  # Non-ground codes
        "bb5": "fst_dens",  # First return density
        "bb8": "1_percentile",
        "bb9": "5_percentile",
        "bba": "25_percentile",
        "bbb": "50_percentile",
        "bbc": "75_percentile",
        "bbd": "95_percentile",
        "bbe": "99_percentile",
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
        #"???": 999,
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
        "bbn": 2024
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
        dict: A dictionary where keys are product codes (e.g., 'bb0') and values are resolution strings.
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


def get_outfnames(
    outputBasename,
    psize=0.5,
    ptile_s=5,
    fpc_psize=10,
    chm_psize=0.2,
    pptiles=(1, 5, 25, 50, 75, 95, 99),
):
    """
    Generate filenames for individually processed tiles and intermediate products.

    Parameters:
        outputBasename (str): Base name for the output files.
        psize (float): Resolution for DEM, maxH, and other products in metres (default: 0.5).
        ptile_s (float): Resolution for percentile tiles in metres (default: 5).
        fpc_psize (float): Resolution for FPC in metres (default: 10).
        chm_psize (float): Resolution for CHM in metres (default: 0.2).
        pptiles (list[int]): List of percentiles to generate filenames for (default: (1, 5, 25, 50, 75, 95, 99)).

    Returns:
        dict: A dictionary where keys are product names and values are their corresponding filenames.
    """
    # Convert resolutions to appropriate string formats
    psize = f"r{int(psize * 100)}cm" if psize < 10 else f"r{int(psize)}m"
    chm_psize = f"r{int(chm_psize * 100)}cm" if chm_psize < 10 else f"r{int(chm_psize)}m"
    ptile_s = f"r{int(ptile_s * 100)}cm" if ptile_s < 10 else f"r{int(ptile_s)}m"
    fpc_psize = f"r{int(fpc_psize * 100)}cm" if fpc_psize < 10 else f"r{int(fpc_psize)}m"

    # Get product stage codes and reverse the dictionary for lookup
    productDict = get_stageDict()
    productNameToCode = {v: k for k, v in productDict.items()}

    fnames = {}

    # List of product names for which to generate filenames
    products_with_psize = [
        'grdR',
        'intens',
        'NonGrd_codes',
        'dem',
        'demHS',
        'maxH',
        'csm',
        'fst_dens',
    ]

    # Generate filenames for products with psize
    for product_name in products_with_psize:
        processing_code = productNameToCode.get(product_name)
        if processing_code is None:
            print(f"Processing code not found for product '{product_name}'")
            continue

        fnames[product_name] = f"{outputBasename}_{processing_code}_{product_name}_{psize}.tif"

    # Handle percentile tiles
    fnames['ptiles'] = {}
    for pos, pptile in enumerate(pptiles):
        product_name = f"{pptile}_percentile"
        processing_code = productNameToCode.get(product_name)
        if processing_code is None:
            print(f"Processing code not found for product '{product_name}'")
            continue

        fnames['ptiles'][pos] = f"{outputBasename}_{processing_code}_{product_name}_{psize}.tif"

    # Handle FPC product
    product_name = 'fpc'
    processing_code = productNameToCode.get(product_name)
    if processing_code is None:
        print(f"Processing code not found for product '{product_name}'")
    else:
        fnames['fpc'] = f"{outputBasename}_{processing_code}_{product_name}_{fpc_psize}.tif"

    # Handle CHM product
    product_name = 'chm'
    processing_code = productNameToCode.get(product_name)
    if processing_code is None:
        print(f"Processing code not found for product '{product_name}'")
    else:
        fnames['chm'] = f"{outputBasename}_{processing_code}_{product_name}_{chm_psize}.tif"

    # Return the filenames dictionary
    return fnames


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
        msg = (
              f"Input {tile} error in filename components.")
        logger.error(msg)
        raise ValueError(msg)
        
    
    # Extract components from the filename
    what = components[0]
    where = components[1]
    when = components[2]
    processing = components[3]

    # Validate the 'what' component
    if what[:2] != "ap":
        msg = (
            f"Filename sensor check failed."
            f"Input {tile} error in filename components.")
        logger.error(msg)
        raise ValueError(msg)

    # Parse the project name
    project = components[4].split(".")[0]
    if len(project[1:]) != 6:
        msg = (
            f"Invalid project name '{project}'. "
            f"Expected 1 prefix character followed by 6 characters.")
        logger.error(msg)
        raise ValueError(msg)
    
    # Create the tile dictionary
    tileDict = {
        "rawtilename": os.path.basename(tile),  # Original tile filename
        "satellite": what[:2],                  # Satellite/platform code (e.g., 'ap')
        "instrument": what[2:4],                # Instrument code
        "returntype": what[4:6],                # Return type (e.g., 'dr')
        "date": int(when),                      # Year of capture
        "tile": where,                          # Tile location (e.g., x448750ys7133000)
        "project": project[1:],                 # Project code (6 characters)
        "stage": processing[:3],                # Processing stage code
        "zone_prefix": processing[-2:-1],       # Zone prefix
        "zoneCode": int(processing[-1:]),       # Zone code
        "xst": np.int32(where[1:7]),            # Starting x-coordinate
        "yst": np.int32(where[9:16]),           # Starting y-coordinate
        "tile_s": np.int32(tile_s),             # Tile size
        "components": components,               # Full list of filename components
    }

    return tileDict



