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
    

########################################################################################################################################
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


########################################################################################################################################

def pdal_convert_epsg():
    """
    just a placeholder for now
    """
    a=a

########################################################################################################################################

def check_las_version(fn,outf,point_format_id=6):
    """
    check las verion - if not version 1.4 then convert
    """
    las_data = laspy.read(fn)
    new_las = laspy.convert(las_data, point_format_id=point_format_id, file_version="1.4")
    new_las.write(outf)
