#!/usr/bin/env python
"""
Display some simple summary statistics for the given project
"""
import argparse
import glob

import numpy as np
import laspy


def getCmdargs():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", help="Directory full of .las/.laz files")
    return p.parse_args()


def main():
    cmdargs = getCmdargs()
    if cmdargs.indir is not None:
        dirname = cmdargs.indir
        filelist = glob.glob(f"{dirname}/*.la[sz]")
    elif cmdargs.infile is not None:
        filelist = [cmdargs.infile]
    else:
        filelist = []

    totalArea = 0
    totalPoints = 0
    densityList = []
    minList = []
    maxList = []
    sizeList = []
    for filename in filelist:
        f = laspy.open(filename)
        hdr = f.header
        xLen = hdr.x_max - hdr.x_min
        yLen = hdr.y_max - hdr.y_min
        area = xLen * yLen
        if area == 0:
            print(f"Warning: file {filename} has zero area")
        nPoints = hdr.point_count
        density = nPoints / area

        totalArea += area
        totalPoints += nPoints
        densityList.append(density)
        minList.append((hdr.x_min, hdr.y_min, hdr.z_min))
        maxList.append((hdr.x_max, hdr.y_max, hdr.z_max))
        sizeList.append((xLen, yLen))

    allMin = np.array(minList).min(axis=0)
    allMax = np.array(maxList).max(axis=0)
    allSize = np.array(sizeList)
    sizeXY = np.median(allSize, axis=0)
    avgPointDensity = totalPoints / totalArea
    densityArr = np.array(densityList)

    crs = hdr.parse_crs()
    epsg = None
    if crs is not None:
        epsg = crs.to_epsg()

    warningMsg = """
        WARNING: Raw LAS/LAZ files are sometimes supplied with errors and omissions in
        their header information, which should be treated with due caution
    """
    print(warningMsg)
    print(f"Projection EPSG: {epsg}")
    print(f"Num files: {len(filelist)}")
    print(f"Total area: {int(round(totalArea / 1E6))} sq-km")
    print(f"Average point density: {avgPointDensity:.1f} pts/sq-m")
    print("Point density distribution (over tiles)")
    print("    Min: {:.1f}".format(densityArr.min()))
    print("    Percentile 25: {:.1f}".format(np.percentile(densityArr, 25)))
    print("    Median: {:.1f}".format(np.median(densityArr)))
    print("    Percentile 75: {:.1f}".format(np.percentile(densityArr, 75)))
    print("    Max: {:.1f}".format(densityArr.max()))
    print("Data ranges")
    print(f"    X: [{allMin[0]:.2f}, {allMax[0]:.2f}]")
    print(f"    Y: [{allMin[1]:.2f}, {allMax[1]:.2f}]")
    print(f"    Z: [{allMin[2]:.2f}, {allMax[2]:.2f}]")
    estimatedTileSize = int(round(sizeXY.mean()))
    print(f"Estimated tile size: {estimatedTileSize}m")


if __name__ == "__main__":
    main()
