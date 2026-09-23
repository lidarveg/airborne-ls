#!/usr/bin/env python
"""
Use RichDEM to calculate flow accumulation
"""
import sys
import os
import argparse

import numpy as np
from osgeo import gdal, gdal_array
from airborne_ls import rw_image_methods

try:
    import richdem
except ImportError:
    cmdname = os.path.basename(sys.argv[0])
    print(f"{cmdname} requires the RichDEM package")
    sys.exit(1)


def getCmdargs():
    p = argparse.ArgumentParser(description="""
            Use RichDEM package to compute stream flow accumulation, with option to
            degrade imagery to lower resolution before computation.
        """)
    p.add_argument("infile", help="Input DEM raster")
    p.add_argument("outfile", help="Output flow accumulation raster")
    p.add_argument("--degradefactor", type=int, default=1,
        help=("Degrade factor (default=%(default)s). " +
              "Input pixel size will be increased by this factor, using averaging"))
    p.add_argument("--method", default="D8",
        choices=['D8', 'D4', 'Rho8', 'Rho4', 'Quinn', 'Freeman', 'Holmgren', 'Dinf'],
        help="Name of method to use, as expected by RichDEM (default=%(default)s)")
    p.add_argument("--driver", default='GTiff',
        help="Output GDAL driver name (default=%(default)s)")
    return p.parse_args()


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    dem = readDEM(cmdargs.infile, cmdargs.degradefactor)
    filledDem = richdem.FillDepressions(dem, epsilon=True)
    flowAcc = richdem.FlowAccumulation(filledDem, method=cmdargs.method).astype(np.float32)
    writeFlowAcc(cmdargs.outfile, flowAcc, cmdargs.driver)


def readDEM(infile, degradefactor):
    """
    Read the input raster, return as a richdem.rdarray
    """
    ds = gdal.Open(infile)
    band = ds.GetRasterBand(1)
    projWKT = ds.GetProjection()
    gt = ds.GetGeoTransform()
    nullval = band.GetNoDataValue()
    (nrows, ncols) = (ds.RasterYSize, ds.RasterXSize)
    # Size of degraded raster
    d = degradefactor
    (nrows_deg, ncols_deg) = (int(nrows // d), int(ncols // d))
    dem = band.ReadAsArray(0, 0, ncols, nrows, ncols_deg, nrows_deg,
        resample_alg=gdal.GRIORA_Average)
    dem_rd = richdem.rdarray(dem, no_data=nullval)
    dem_rd.projection = projWKT
    gt_deg = (gt[0], gt[1] * d, 0, gt[3], 0, gt[5] * d)
    dem_rd.geotransform = gt_deg
    dem_rd.metadata = {}
    dem_rd.metadata.update(ds.GetMetadata())

    return dem_rd


def writeFlowAcc(outfile, flowAcc, driverName):
    """
    Write the flow accumulation to a file
    """
    drvr = gdal.GetDriverByName(driverName)
    if os.path.exists(outfile):
        drvr.Delete(outfile)
    (nrows, ncols) = flowAcc.shape
    gdalType = gdal_array.NumericTypeCodeToGDALTypeCode(flowAcc.dtype)
    creationoptions = rw_image_methods.creationOptionsByDriver.get(driverName, [])
    ds = drvr.Create(outfile, ncols, nrows, 1, eType=gdalType, options=creationoptions)
    band = ds.GetRasterBand(1)
    band.WriteArray(flowAcc)
    band.SetNoDataValue(flowAcc.no_data)
    band.ComputeStatistics(approx_ok=False)
    ds.SetGeoTransform(flowAcc.geotransform)
    ds.SetProjection(flowAcc.projection)
    ds.BuildOverviews(resampling="AVERAGE", overviewlist=[2, 4, 8, 16, 32])


if __name__ == "__main__":
    main()
