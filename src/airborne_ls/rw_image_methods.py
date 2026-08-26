#!/usr/bin/env python

"""
Script for processing LAS/LAZ files using laspy and converting them to a structured
numpy record array.
"""
import numpy as np
from osgeo import gdal, osr, gdal_array


gdal.UseExceptions()


creationOptionsByDriver = {
    "GTiff": ["COMPRESS=DEFLATE", "TILED=YES", "INTERLEAVE=BAND", "BIGTIFF=IF_SAFER"],
    "HFA": ["COMPRESSED=YES", "IGNOREUTM=YES"],
    "KEA": [],
    # COG is only used for the mosaics
    "COG": ["COMPRESS=DEFLATE", "BLOCKSIZE=256", "RESAMPLING=BILINEAR", "BIGTIFF=IF_SAFER"]
}


###################################################################################################


def writeImage(image, outfile, driverName="GTiff", tlx=0.0, tly=0.0,
               binsize=0.0, epsg=None, nullVal=None, parent_file=None,
               overviewResampling="BILINEAR", overviewLevels=[4, 8, 16, 32, 64],
               creationoptions=None):
    """
    Write data to a GDAL supported image file format

    Parameters:
    - image: The image data array.
    - outfile: The output filename.
    - driverName: The GDAL driver to use.
    - tlx, tly: Top-left x and y coordinates.
    - binsize: Pixel size.
    - epsg: EPSG code for spatial reference.
    - nullVal: NoData value.
    - parent_file: The parent LAS file used to generate the image.


    GTiff Driver Issues: When using the GTiff driver, I encountered errors due to
    mismatches between the file extension and the driver or due to incorrect creation options.

    To resolve this issue, we need to switching to the GTiff driver and ensuring that:
    The Output File Extension Matches the Driver: When using the GTiff driver, the output
    file should have a .tif extension, not .img.

    """
    if len(image.shape) == 2:
        ny, nx = image.shape
        nz = 1
    if len(image.shape) == 3:
        nz, ny, nx = image.shape
    drvr = gdal.GetDriverByName(driverName)
    dt = image.dtype

    # Map numpy dtype to GDAL data type
    gdaldtype = gdal_array.NumericTypeCodeToGDALTypeCode(dt)
    if gdaldtype is None:
        raise ValueError(f"Unsupported data type: {dt}")

    if creationoptions is None:
        creationoptions = creationOptionsByDriver.get(driverName, [])

    ds = drvr.Create(outfile, nx, ny, nz, gdaldtype, options=creationoptions)
    ds.SetGeoTransform([tlx, binsize, 0, tly, 0, -binsize])

    if epsg is not None:
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(epsg)
        ds.SetProjection(proj.ExportToWkt())

    if nz > 1:
        for i in range(nz):
            band = ds.GetRasterBand(i + 1)
            band.WriteArray(image[i, :, :], 0, 0)
    else:
        band = ds.GetRasterBand(1)
        band.WriteArray(image, 0, 0)

    # Set the null value on every band
    if nullVal is not None:
        for i in range(nz):
            band = ds.GetRasterBand(i + 1)
            band.SetNoDataValue(nullVal)

    ds.BuildOverviews(resampling=overviewResampling, overviewlist=overviewLevels)

    ds.FlushCache()
    band = None
    ds = None


#############################################################################################
# image rw
def imgH(img):
    """
    read imagine header and return as a dictionary
    """
    data = gdal.Open(img)
    info = np.array(data.GetGeoTransform())
    bands = data.RasterCount
    xst = info[0]
    yst = info[3]
    xdim = data.RasterXSize
    ydim = data.RasterYSize
    pixel_s = info[1]
    # proj = data.GetProjection()
    h = {
        "xdim": xdim,
        "ydim": ydim,
        "bands": bands,
        "pixel_s": pixel_s,
        "tlx": xst,
        "tly": yst,
    }

    return h


def imgRead(infile):
    """
    read image via gdal
    """
    oFile = gdal.Open(infile)
    img = oFile.GetRasterBand(1).ReadAsArray()
    return img


#############################################################################################
def readtxt(fn):
    """
    read comma seperated txt file
    """
    with open(fn) as f:
        txt = [line.strip() for line in f]
    col = len((txt[0]).split(","))
    row = len(txt)
    data = np.zeros((col, row), dtype=np.float64)
    for ct, row in enumerate(txt):
        data[:, ct] = row.split(",")
    return np.array(data)


def write_txt_all(data, outfile):
    """
    write comma seperated txt file
    """
    if ((data.shape)[1]) >= (data.shape)[0]:
        data = np.transpose(data)

    with open(outfile, "w") as fout:
        for val in range(int((data.shape)[0])):
            res = ",".join(str(item) for item in data[val, :])
            line_to_write = f"{res}\n"
            fout.write(line_to_write)
    fout.close()


#############################################################################################


def force_header(header, tile_s):
    """
    Adjust the header boundaries to ensure all arrays are aligned to a grid of size tile_s x tile_s.

    Parameters:
        header (object): The header object containing x_min, x_max, y_min, and y_max attributes.
        tile_s (float): The size of the tile to align the header boundaries to.

    Returns:
        object: The updated header with adjusted boundaries.
    """
    header.x_min = (header.x_min // tile_s) * tile_s
    header.x_max = tile_s + (header.x_max // tile_s) * tile_s
    header.y_min = (header.y_min // tile_s) * tile_s
    header.y_max = np.ceil(header.y_max / tile_s) * tile_s
    return header


def get_mmXYZ(x, y, z):
    """
    Extract the minimum and maximum values for X, Y, and Z coordinates.

    Parameters:
        x (array-like): Array of X coordinates.
        y (array-like): Array of Y coordinates.
        z (array-like): Array of Z coordinates.

    Returns:
        tuple: A tuple containing (minX, maxX, minY, maxY, minZ, maxZ).
    """
    minX = np.min(x)
    maxX = np.max(x)
    minY = np.min(y)
    maxY = np.max(y)
    minZ = np.min(z)
    maxZ = np.max(z)

    return minX, maxX, minY, maxY, minZ, maxZ


###################################################################################################


def setColorTable(imgfile, clrTblArr):
    """
    Set the colour table on the given Byte image file.

    The colour table is given as a 2D array of shape (256, 3). Each row is
    a colour. The i-th row is the red/green/blue values for the pixel value i.
    The RGB values are integers in the range [0, 255].

    The colour table is set only on the first band of the image file.

    Parameters:
      imgfile (str): Name of image file
      clrTblArr (2-d array): Colour table values
    """
    if clrTblArr.shape != (256, 3):
        raise ValueError("clrTblArr must be shape (256, 3)")

    ds = gdal.Open(imgfile, gdal.GA_Update)
    bandobj = ds.GetRasterBand(1)

    clrTbl = gdal.ColorTable()
    for i in range(len(clrTblArr)):
        colEntry = tuple(clrTblArr[i])
        clrTbl.SetColorEntry(i, colEntry)
    bandobj.SetRasterColorTable(clrTbl)
    bandobj.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
