"""
Generate a test file of Lidar point data. This can then be used to run unit tests
of the various components of the package.

"""
import numpy
import laspy
import pyproj

from airborne_ls import const


def generateTestPointData(filename):
    """
    Generate a single file of Lidar point data.

    Parameters:
      filename (str): Name of file in which to write generated data
    """
    epsg = 28356    # UTM56, AMG

    # A regular grid of (x, y) points. The values below give us 4 points / square metre
    n = 2000
    step = 0.5
    (xMin, yMin) = (485000 + step / 2, 6963000 + step / 2)
    (xMax, yMax) = (xMin + n * step, yMin + n * step)
    # The (x, y) coordinates of every point
    (y, x) = numpy.mgrid[yMin:yMax:step, xMin:xMax:step]
    y = numpy.flip(y, axis=0)

    # Assign some Z values to each point.
    # Begin with flat ground at 100m, with a bilinear smooth fall of 2m, from the
    # bottom-left to the top-right. We want a small fall so our flow accumulation
    # will not be completely flummoxed.
    drop = 2
    hgtBL = 100
    hgtTL = hgtBL - drop / 2
    hgtBR = hgtBL - drop / 2
    hgtTR = hgtBL - drop
    # Bilinear interpolation between these 4 corners
    z = (
        hgtBL * ((xMax - x) / (xMax - xMin)) * ((yMax - y) / (yMax - yMin)) +
        hgtBR * ((x - xMin) / (xMax - xMin)) * ((yMax - y) / (yMax - yMin)) +
        hgtTL * ((xMax - x) / (xMax - xMin)) * ((y - yMin) / (yMax - yMin)) +
        hgtTR * ((x - xMin) / (xMax - xMin)) * ((y - yMin) / (yMax - yMin))
    )

    # Add a conical hill
    (hillCtrX, hillCtrY) = (485280, 6963750)
    hillRadius = 75
    hMax = 50
    d = numpy.sqrt((x - hillCtrX)**2 + (y - hillCtrY)**2)
    hgt = hMax * (1 - d / hillRadius)
    hgtGt0 = (hgt > 0)
    z[hgtGt0] = z[hgtGt0] + hgt[hgtGt0]

    classification = numpy.full(x.shape, const.PTCLASS_GROUND, dtype=numpy.uint8)

    # Add a gully. It starts at the bottom of the hill (with Z == 100m) and runs east, down
    # to z == (100 - endDepth)
    (gullyStartX, gullyStartY) = (hillCtrX + hillRadius, hillCtrY)
    (gullyEndX, gullyEndY) = (xMax, gullyStartY)
    gullyLen = gullyEndX - gullyStartX
    endDepth = 30
    # Slope of gully sides (running in Y direction)
    slopeY = numpy.radians(15)
    tanSlopeY = numpy.tan(slopeY)
    slopeRun = endDepth / tanSlopeY
    # A mask for the rectangle surrounding the gully
    gullySthEdge = gullyEndY - slopeRun
    gullyNthEdge = gullyEndY + slopeRun
    gullyRectangleMask = ((x > gullyStartX) & (x <= gullyEndX) &
                         (y > gullySthEdge) & (y <= gullyNthEdge))
    # Masks for the north and south slopes of the gully
    gullySthMask = (gullyRectangleMask & (y < gullyEndY))
    gullyNthMask = (gullyRectangleMask & (y >= gullyEndY))

    # Compute depth for each slope
    gullyDepthSth = ((endDepth * (y - gullySthEdge) / slopeRun) *
                     (x - gullyStartX) / gullyLen).clip(0, endDepth)
    gullyDepthNth = ((endDepth * (gullyNthEdge - y) / slopeRun) *
                     (x - gullyStartX) / gullyLen).clip(0, endDepth)
    # Subtract from Z
    z[gullySthMask] = z[gullySthMask] - gullyDepthSth[gullySthMask]
    z[gullyNthMask] = z[gullyNthMask] - gullyDepthNth[gullyNthMask]

    # Add a circular lake in the south-west quadrant. This will have sloping sides and a
    # flat bottom, and all bottom points (i.e. lake surface) will be classified as water
    # returns instead of ground.
    (lakeCtrX, lakeCtrY) = (485280, 6963250)
    lakeRadius = 100
    depthMax = 50
    d = numpy.sqrt((x - lakeCtrX)**2 + (y - lakeCtrY)**2)
    depth = depthMax * (1 - d / lakeRadius)
    depthGt0 = (depth > 0)
    middleOfLake = (d < 0.75 * lakeRadius)
    depthToWater = depth[middleOfLake].min()
    depth[middleOfLake] = depthToWater
    z[depthGt0] = z[depthGt0] - depth[depthGt0]
    # Convert points in the middle to water returns
    classification[middleOfLake] = const.PTCLASS_WATER

    # Create the output file and write the points
    header = laspy.LasHeader(point_format=1, version="1.4")
    header.offsets = numpy.array([0.0, 0.0, 0.0])
    header.scales = numpy.array([0.01, 0.01, 0.01])
    # Include projection
    crsObj = pyproj.CRS.from_epsg(epsg)
    header.add_crs(crsObj)

    pts = laspy.LasData(header)
    pts.x = x.flatten()
    pts.y = y.flatten()
    pts.z = z.flatten()
    pts.classification = classification.flatten()

    pts.write(filename)
