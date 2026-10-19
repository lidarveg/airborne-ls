#!/usr/bin/env python
"""
Re-tile a set of LAS/LAZ files into square tiles of a fixed size.
"""
import os
import argparse
import glob
import copy

import numpy as np
import laspy

from airborne_ls import lazfile_rw


def getCmdargs():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser(description="""
            Re-tile a set of LAS/LAZ files into square tiles of a fixed size. This is most
            commonly useful when the starting point is un-tiled flightline files, because
            the rest of the ALSV processing chain assumes small square tiles of point data.
            The output tiles are written uncompressed,, and intended to be temporary"
        """)
    p.add_argument("--indir", required=True,
        help="Directory containing LAS/LAZ files to process.")
    p.add_argument("--outdir", required=True,
        help="Directory for output LAS tiles. Must be different to --indir")
    p.add_argument("--tilesize", type=int, default=1000,
        help="Size of output tiles in metres (default=%(default)s)")
    p.add_argument("--noparallelread", default=False, action="store_true",
        help=("Disable laspy's parallel reading (& writing) of LAS/LAZ files. By default, " +
              "laspy will try to use as many CPU cores as it can"))
    p.add_argument("--minlasversion", default="1.4",
        help=("Minimum LAS format version for output files (default=%(default)s). " +
              "Newer format input files will preserve their newer version"))

    cmdargs = p.parse_args()

    if cmdargs.outdir == cmdargs.indir:
        msg = "--outdir must be different to --indir"
        raise ValueError(msg)

    if not os.path.exists(cmdargs.outdir):
        msg = f"Out directory '{cmdargs.outdir}' does not exist"
        raise ValueError(msg)
    if not os.path.isdir(cmdargs.outdir):
        msg = f"Out directory '{cmdargs.outdir}' is not a directory"
        raise ValueError(msg)

    lasVersNum = tuple([int(i) for i in cmdargs.minlasversion.split('.')])
    if len(lasVersNum) > 2:
        msg = f"LAS version {cmdargs.minlasversion} has too many components ({len(lasVersNum)})"
        raise ValueError(msg)
    cmdargs.minlasversion = laspy.header.Version(*lasVersNum)

    return cmdargs


def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()

    laspyBackend = lazfile_rw.chooseLaspyBackend(cmdargs.noparallelread)

    fileList = sorted(glob.glob(f"{cmdargs.indir}/*.la[sz]"))
    for filename in fileList:
        print("Retiling", filename)
        data = laspy.read(filename, laz_backend=laspyBackend)
        # Convert older formats to the requested LAS version.
        if data.header.version < cmdargs.minlasversion:
            data = laspy.convert(data, file_version=cmdargs.minlasversion)

        pointsTileX = np.int64(np.floor(data.x / cmdargs.tilesize) * cmdargs.tilesize)
        pointsTileY = np.int64(np.ceil(data.y / cmdargs.tilesize) * cmdargs.tilesize)

        possibleXtiles = np.unique(pointsTileX)
        possibleYtiles = np.unique(pointsTileY)

        for tileY in possibleYtiles:
            for tileX in possibleXtiles:
                ndxInTile = ((tileY == pointsTileY) & (tileX == pointsTileX))
                numPoints = np.count_nonzero(ndxInTile)
                if numPoints > 0:
                    pointsInTile = data[ndxInTile]
                    outFile = f"{cmdargs.outdir}/retiled_x{tileX}y{tileY}.las"

                    # Cope with any left-over empty files. Mainly coping with error states.
                    if os.path.exists(outFile):
                        fileSize = os.stat(outFile).st_size
                        if fileSize == 0:
                            os.remove(outFile)

                    if not os.path.exists(outFile):
                        new_hdr = copy.copy(data.header)
                        new_hdr.point_count = 0
                        new_las = laspy.LasData(new_hdr)
                        new_las.points = pointsInTile.points.copy()
                        new_las.write(outFile)
                        del new_las
                    else:
                        with laspy.open(outFile, mode='a') as outlas:
                            outlas.append_points(pointsInTile.points.copy())
                    del pointsInTile
                del ndxInTile


if __name__ == "__main__":
    main()
