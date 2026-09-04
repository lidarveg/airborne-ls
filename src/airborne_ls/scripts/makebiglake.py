#!/usr/bin/env python
"""
For testing purposes only. This script takes two adjacent DEM tiles, and imposes
a null-filled hole crossing the boundary between the two. This is intended to
represent a water body where the bin-level interpolation failed, leaving this
substantial hole. This can be used to test whether the gap-filling script is
able to interpolate across such a hole.

Given that it is intended to represent a water body, this script makes a substantial
effort to create a hole with its edges all at the same elevation (a characteristic
property of a water body). In pursuit of that goal, it also creates a shoreline
which slopes upwards from that level to join the genuine terrain at some distance.

"""
import argparse

import numpy
from osgeo import gdal


gdal.UseExceptions()


def getCmdargs():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser(description="""
        Carve an artificial hole into terrain across two adjacent tiles, representing
        a large water body. The area corresponding to water is filled with nulls. The two
        tiles must be adjacent in the east-west direction.
    """)
    p.add_argument("westtile", help="Western tile")
    p.add_argument("easttile", help="Eastern tile")
    p.add_argument("--nullval", type=float, default=-999.9,
        help="Null value to fill water area (default=%(default)s)")
    cmdargs = p.parse_args()
    return cmdargs


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    (westtile, nullVal) = readDem(cmdargs.westtile)
    (easttile, _) = readDem(cmdargs.easttile)

    dem = numpy.concatenate([westtile, easttile], axis=1)
    (nRows, nCols) = dem.shape
    tileNcols = nCols // 2

    # The dimensions of the inner and outer rectangles, expressed in the row and column
    # numbers of the combined array. The inner rectangle is the lake surface, the outer
    # rectangle is the top of the shore slope to the genuine terrain.
    slopeWidth = 100
    (iLen, oLen) = (tileNcols - slopeWidth, tileNcols + slopeWidth)
    (iWidth, oWidth) = (nRows // 10 - slopeWidth, nRows // 10 + slopeWidth)
    (oLeft, oTop) = (tileNcols // 2, nRows // 2 - oWidth // 2)
    (iLeft, iTop) = (oLeft + slopeWidth, oTop + slopeWidth)
    (oRight, oBottom) = (oLeft + oLen, oTop + oWidth)
    (iRight, iBottom) = (iTop + iLen, iTop + iWidth)

    # Extract pixels on outer boundary
    outerPixels = numpy.concatenate([
        dem[oTop, oLeft:oRight], dem[oBottom, oLeft:oRight],
        dem[oTop + 1:oBottom - 1, oLeft], dem[oTop + 1:oBottom - 1, oRight]
    ])
    # Extract pixels in the slope region
    slopePixels = numpy.concatenate([
        dem[oTop:iTop, oLeft:oRight].flatten(),
        dem[iBottom:oBottom, oLeft:oRight].flatten(),
        dem[iTop + 1:iBottom - 1, oLeft:iLeft].flatten(),
        dem[iTop + 1:iBottom - 1, iRight:oRight].flatten()
    ])
    # Some extreme elevation values from the boundary and slope area
    boundaryMax = outerPixels.max()
    slopeMin = slopePixels.min()
    # The elevation step size to reduce for a one-pixel step
    zStep = (boundaryMax - slopeMin) / slopeWidth

    # Now excavate the slope
    for i in range(slopeWidth):
        z = boundaryMax - i * zStep
        # Define an intermediate rectangle, <i> pixels in from the outer one
        top = oTop + i
        bottom = oBottom - i
        left = oLeft + i
        right = oRight - i
        # On this rectangle, any pixels greater than z are set to z
        topRowSlice = (top, slice(left, right))
        bottomRowSlice = (bottom, slice(left, right))
        leftColSlice = (slice(top, bottom + 1), left)
        rightColSlice = (slice(top, bottom + 1), right)
        for s in [topRowSlice, bottomRowSlice, leftColSlice, rightColSlice]:
            # Locations within dem slice which are greater than z
            ndxHigh = dem[s] > z
            # Set these to z
            dem[s][ndxHigh] = z
    # Inside the inner rectangle, fill with nulls
    dem[top:bottom, left:right] = nullVal

    rewriteDem(cmdargs.westtile, dem[:, :tileNcols])
    rewriteDem(cmdargs.easttile, dem[:, tileNcols:])


def readDem(filename):
    """
    Read the DEM tile into an array. Return the array, and the file's null value
    """
    ds = gdal.Open(filename)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    nullVal = band.GetNoDataValue()
    return (arr, nullVal)


def rewriteDem(filename, dem):
    """
    Rewrite the given DEM array into the given file
    """
    ds = gdal.Open(filename, gdal.GA_Update)
    band = ds.GetRasterBand(1)
    band.WriteArray(dem)
    ds.BuildOverviews(resampling="BILINEAR", overviewlist=[4, 8, 16, 32, 64])


if __name__ == "__main__":
    main()
