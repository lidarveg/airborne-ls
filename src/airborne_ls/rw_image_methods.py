#!/usr/bin/env python

"""
Script for processing LAS/LAZ files using laspy and converting them to a structured
numpy record array.
"""
import numpy as np
from osgeo import gdal, osr, gdal_array

# from airborne-ls import


gdal.UseExceptions()


###################################################################################################


def writeImage(image, outfile, cmdargs, driver="GTiff", tlx=0.0, tly=0.0,
               binsize=0.0, epsg=None, nullVal=None, parent_file=None,
               overviewResampling="BILINEAR", overviewLevels=[4, 8, 16, 32, 64]):
    """
    Write data to a GDAL supported image file format

    Parameters:
    - image: The image data array.
    - outfile: The output filename.
    - driver: The GDAL driver to use.
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
    driver = gdal.GetDriverByName(driver)
    dt = image.dtype

    # Map numpy dtype to GDAL data type
    gdaldtype = gdal_array.NumericTypeCodeToGDALTypeCode(dt)
    if gdaldtype is None:
        raise ValueError(f"Unsupported data type: {dt}")

    ds = driver.Create(outfile, nx, ny, nz, gdaldtype)
    ds.SetGeoTransform([tlx, binsize, 0, tly, 0, -binsize])

    if epsg is not None:
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(epsg)
        ds.SetProjection(proj.ExportToWkt())

    # set colour table for bb4
    if outfile.find("NonGrdCodes") > 0:
        colors = gdal.ColorTable()
        # set color for each value
        colors.SetColorEntry(0, (254, 254, 254))  # never classified:
        colors.SetColorEntry(1, (200, 200, 200))  # unclassified: light gray
        colors.SetColorEntry(2, (0, 0, 0))        # ground classification
        colors.SetColorEntry(3, (0, 240, 0))      # low veg: green1
        colors.SetColorEntry(4, (0, 160, 0))      # medium veg: green2
        colors.SetColorEntry(5, (0, 80, 0))       # high veg: green3
        colors.SetColorEntry(6, (255, 0, 0))      # building: red
        colors.SetColorEntry(7, (255, 255, 0))
        colors.SetColorEntry(8, (255, 255, 0))
        colors.SetColorEntry(9, (0, 0, 255))      # blue for water
        colors.SetColorEntry(10, (255, 0, 255))
        colors.SetColorEntry(11, (255, 20, 255))
        colors.SetColorEntry(12, (255, 30, 255))
        colors.SetColorEntry(13, (255, 40, 255))
        colors.SetColorEntry(14, (255, 50, 255))
        colors.SetColorEntry(15, (255, 60, 255))
        colors.SetColorEntry(16, (255, 70, 255))
        colors.SetColorEntry(17, (255, 80, 255))
        colors.SetColorEntry(18, (255, 90, 255))
        colors.SetColorEntry(19, (255, 100, 255))
        colors.SetColorEntry(254, (101, 67, 33))  # brown for background

    if nz > 1:
        for i in range(nz):
            band = ds.GetRasterBand(i + 1)
            # set color table and color interpretation
            if outfile.find("NonGrdCodes") > 0:
                band.SetRasterColorTable(colors)
                band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
            band.WriteArray(image[i, :, :], 0, 0)
    else:
        band = ds.GetRasterBand(1)
        # set color table and color interpretation
        if outfile.find("NonGrdCodes") > 0:
            band.SetRasterColorTable(colors)
            band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
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
