#!/usr/bin/env python
"""
Run tests of major functionality, using generated input data
"""
import sys
import os
import argparse
import subprocess
import tempfile
import shutil
import glob

import numpy as np
from osgeo import gdal
import laspy
from airborne_ls.tests import gentestdata
from airborne_ls import qvf


gdal.UseExceptions()


def getCmdargs():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser()
    p.add_argument("--keep", default=False, action="store_true",
        help="Keep test files after completion (default will delete them all)")
    return p.parse_args()


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    tmpdir = tempfile.mkdtemp(prefix='alsv_tests_', suffix='.tmpd', dir='.')
    lazfile = os.path.join(tmpdir, 'testpoints.laz')
    indexedDir = os.path.join(tmpdir, 'indexedlaz')

    testCounts = TestCounts()
    createTestData(lazfile, testCounts)
    checkStandardisation(lazfile, indexedDir, testCounts)
    checkTileProd(indexedDir, testCounts)
    checkGapFill(indexedDir, testCounts)
    checkFlowAcc(indexedDir, testCounts)

    if not cmdargs.keep:
        shutil.rmtree(tmpdir)

    if testCounts.failCount > 0:
        print(f"Failed {testCounts.failCount} tests, passed {testCounts.passCount}")
        sys.exit(1)
    else:
        print(f"Passed all {testCounts.passCount} tests")


def createTestData(lazfile, testCounts):
    """
    Generate test data and check it
    """
    testName = "GeneratePoints"
    gentestdata.generateTestPointData(lazfile)

    f = laspy.open(lazfile)
    hdr = f.header

    checkEqual(testName, hdr.point_count, 4015785, testCounts, "Generated point count")

    crs = hdr.parse_crs()
    if crs is None:
        reportError(testName, "No CRS in test file")
        testCounts.failed()
    else:
        testCounts.passed()

    if crs is not None:
        epsg = crs.to_epsg()
        checkEqual(testName, epsg, 28356, testCounts, "EPSG")


def checkStandardisation(lazfile, indexedDir, testCounts):
    """
    Test standardisation of LAZ files
    """
    testName = "StandardiseLidar"
    if not os.path.exists(lazfile):
        reportError(testName, f"Skipped, {lazfile} not found")
        testCounts.failed()

    if not os.path.exists(indexedDir):
        os.mkdir(indexedDir)

    projectName = "testpr"
    cmd = ['alsv_lidar_standardisation', '--infile', lazfile, '--outdir', indexedDir,
           '--intilesize', '1000', '--ii', 'uk', '--project', projectName, '--year', '2025']
    ok = runCmd(cmd, testName)

    if ok:
        stdLazfile = None
        outfileList = glob.glob(f"{indexedDir}/*.laz")
        if len(outfileList) == 0:
            reportError(testName, "No standardised LAZ file found")
            testCounts.failed()
        elif len(outfileList) > 1:
            reportError(testName, f"Found {len(outfileList)} output LAZ files")
            testCounts.failed()
        else:
            stdLazfile = outfileList[0]

        if stdLazfile is not None:
            f = laspy.open(lazfile)
            inHdr = f.header
            f = laspy.open(outfileList[0])
            stdHdr = f.header

            checkEqual(testName, stdHdr.point_count, inHdr.point_count, testCounts,
                "Output point count")
    else:
        testCounts.failed()


def checkTileProd(indexedDir, testCounts):
    """
    Check the raster tile products
    """
    testName = "TileProducts"

    cmd = ['alsv_tile_products', '--indir', indexedDir, '--tilesize', '1000']
    ok = runCmd(cmd, testName)

    if ok:
        checkDEM(indexedDir, testCounts)
        checkFPC(indexedDir, testCounts)
    else:
        testCounts.failed()


def checkDEM(indexedDir, testCounts):
    """
    Find the DEM file and do some basic checks on it
    """
    testName = "checkDEM"
    nullVal = -999.0
    demfileList = glob.glob(f"{indexedDir}/*/*_bb0m6_ldem_*.tif")
    if len(demfileList) == 0:
        reportError(testName, "DEM file not found")
        testCounts.failed()
    else:
        dem = readImg(demfileList[0])
        nullCount = np.count_nonzero(dem == nullVal)
        (minHgt, maxHgt) = (dem[dem != nullVal].min(), dem.max())
        checkEqual(testName, minHgt, 68.34, testCounts, "Min hgt")
        checkEqual(testName, maxHgt, 148.74, testCounts, "Max hgt")
        checkEqual(testName, nullCount, 50338, testCounts, "Null count")


def checkFPC(indexedDir, testCounts):
    """
    Find the FPC image and do some simple checks
    """
    testName = "checkFPC"
    fpcfileList = glob.glob(f"{indexedDir}/*/*_bbhm6_lfpc_*.tif")
    if len(fpcfileList) == 0:
        reportError(testName, "FPC file not found")
        testCounts.failed()
    else:
        fpc = readImg(fpcfileList[0])
        numZero = np.count_nonzero(fpc == 0)
        nonzeroFPC = fpc[fpc > 0]
        forestMeanFPC = nonzeroFPC.mean()
        checkEqual(testName, numZero, 9952, testCounts, "Number of nonzero FPC pixels")
        checkEqual(testName, forestMeanFPC, 43.333333333333336, testCounts,
            "Forest FPC mean")


def checkGapFill(indexedDir, testCounts):
    """
    Check the dem_gapfill command
    """
    testName = "GapFillDEM"
    nullVal = -999.0

    cmd = ['alsv_dem_gapfill', '--indir', indexedDir]
    ok = runCmd(cmd, testName)
    if ok:
        gfdemList = glob.glob(f"{indexedDir}/*/*_bb0m6_lgapfilleddem_*.tif")
        if len(gfdemList) == 0:
            reportError(testName, "No gap-filled dem found")
            testCounts.failed()
        else:
            dem = readImg(gfdemList[0])
            nullCount = np.count_nonzero(dem == nullVal)
            (minHgt, maxHgt) = (dem[dem != nullVal].min(), dem.max())
            checkEqual(testName, minHgt, 68.34, testCounts, "Min hgt")
            checkEqual(testName, maxHgt, 148.74, testCounts, "Max hgt")
            checkEqual(testName, nullCount, 0, testCounts, "Null count")
    else:
        testCounts.failed()


def checkFlowAcc(indexedDir, testCounts):
    """
    Check flow accumulation command
    """
    testName = "FlowAcc"

    gfdemList = glob.glob(f"{indexedDir}/*/*_bb0m6_lgapfilleddem_*.tif")
    if len(gfdemList) == 1:
        gfdemfile = gfdemList[0]
    else:
        gfdemfile = None

    if gfdemfile is not None and os.path.exists(gfdemfile):
        flowaccfile = qvf.setoptionfield(gfdemfile, 'l', 'flowacc')
        cmd = ['alsv_flowacc', gfdemfile, '--outfile', flowaccfile]
        ok = runCmd(cmd, testName)
        if ok:
            if os.path.exists(flowaccfile):
                flowacc = readImg(flowaccfile)
                # Rows and columns of pixels with largest flow
                (rows, cols) = np.where(flowacc > 1000000)
                uniqRows = np.unique(rows)
                checkEqual(testName, len(uniqRows), 1, testCounts, "Large flow row count")
                checkEqual(testName, uniqRows[0], 499, testCounts, "Large flow row")
                checkEqual(testName, len(cols), 675, testCounts, "Large flow length")
            else:
                reportError(testName, "Flow accumulation file not found")
                testCounts.failed()
        else:
            testCounts.failed()
    else:
        reportError(testName, "Gap-filled DEM file '{gfdemfile}' not found")
        testCounts.failed()


def runCmd(cmd, testName):
    """
    Run a command in a subprocess, and return True if all OK.

    If error detected, report with reportError().

    Parameters:
      cmd (list of str): List of command and args, suitable for Popen

    Returns:
      ok (bool): True if return code == 0 and len(stderr) == 0
    """
    ok = True
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True)
        (stdout, stderr) = proc.communicate()
        exitStat = proc.returncode
        # Remove nasty cursor control string that RichDEM puts into stderr
        badStr = '\x1b[2K'      # '<esc>[2K', erase line
        stderr = stderr.replace(badStr, '').strip()

        if exitStat != 0 or len(stderr) > 0:
            msg = f"Exit status {exitStat} from {cmd[0]}.\n{stderr}"
            reportError(testName, msg)
            ok = False
    except FileNotFoundError:
        reportError(testName, f"Command '{cmd[0]}' not found")
        ok = False
    except Exception as e:
        reportError(testName, str(e))
        ok = False

    return ok


def readImg(filename):
    """
    Read band 1 of given image file, return as 2-d array
    """
    ds = gdal.Open(filename)
    img = ds.GetRasterBand(1).ReadAsArray()
    return img


def checkEqual(testName, val1, val2, testCounts, msg):
    """
    Check the two values are equal, report if not. Update testCounts.
    """
    if val1 != val2:
        reportError(testName, f"{msg} {val1} != {val2}")
        testCounts.failed()
    else:
        testCounts.passed()


def reportError(testName, msg):
    print(f"{testName}: {msg}")


class TestCounts:
    def __init__(self):
        self.passCount = 0
        self.failCount = 0

    def passed(self):
        self.passCount += 1

    def failed(self):
        self.failCount += 1


if __name__ == "__main__":
    main()
