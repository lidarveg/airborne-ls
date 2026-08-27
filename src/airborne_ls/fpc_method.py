#!/usr/bin/env python

"""
coding of method:

    Fisher, A., Armston, J., Goodwin, N., Scarth, P. (2020). Modelling canopy gap probability,
    foliage projective cover and crown projective cover from airborne lidar metrics in
    Australian forests and woodlands. Remote Sensing of Environment, 237, 111520.
    doi: 10.1016/j.rse.2019.111520

"""

import numpy as np
from numba import njit

from airborne_ls import gridding_methods


def fix_pulse_positions(data):
    """
    This changes the x and y coordinates of subsequent returns to that of the
    first return, to make sure they are in the same pixel.
    """
    nElems = len(data)
    ret_diff = data["RETURN_NUMBER"][1:nElems] - data["RETURN_NUMBER"][0 : nElems - 1]
    locs = np.argwhere(ret_diff == 1)
    data["X"][locs + 1] = data["X"][locs]
    data["Y"][locs + 1] = data["Y"][locs]

    return data


def check_pts_pulses(las_data):
    """
    Assigns flight line numbers to points based on timestamp differences.
    """
    nElems = len(las_data)
    flightLines = np.zeros(nElems, dtype=np.uint8)
    diff = las_data["TIMESTAMP"][1:] - las_data["TIMESTAMP"][:-1]
    idx = np.argwhere(diff > 10) + 1

    if len(idx) > 0:
        st = 0
        for ct, v in enumerate(idx):
            flightLines[st : v[0]] = ct
            st = v[0]
        flightLines[st:] = len(idx) + 1

    return flightLines


@njit
def fpcGridding(row, col, hgt, wgt, fpc, canopyThreshold):
    """
    Calculate two bands:
     - sum of return weights above canopy thershold (canopy weights)
     - sum of all return weights (total weights)
    """
    w2 = np.copy(wgt)
    w2[hgt < canopyThreshold] = 0.0
    for p in range(len(hgt)):
        fpc[0, row[p], col[p]] += w2[p]
        fpc[1, row[p], col[p]] += wgt[p]


def doFPC(xMin, yMax, data, flightLines, heightAboveGround, fpc_size, tile_s,
        split_fpc, canopyThreshold):
    """
    Calculates FPC as the proportion of weighted returns above the canopy threshold.
    FPC is calculated separately for different flight lines and combined using the mean.
    All returns from each pulse are located in the pixel of the first return.

    Some datasets have very high point densities at the edge of scans as the mirror
    changes direction and slows down. No perfect solution as yet.
    """
    data = fix_pulse_positions(data)
    nRows = int(np.ceil(tile_s / fpc_size))
    nCols = nRows  # Assumes square tiles.
    nullVal = 255

    if split_fpc:
        flights = np.unique(flightLines)
        numFlights = flights.size
        fpc_flights = np.zeros((numFlights, nRows, nCols), dtype=np.float32)
        fpcFP = np.zeros((2, nRows, nCols), dtype=np.float32)

        for f, flight in enumerate(flights):
            idx = flightLines == flight
            x = np.copy(data["X"][idx])
            y = np.copy(data["Y"][idx])
            hag = np.copy(heightAboveGround[idx])
            numberOfReturns = data["NUMBER_OF_RETURNS"][idx]

            weight = 1 / numberOfReturns.astype(np.float32)
            row, col = gridding_methods.xyToRowCol(x, y, xMin, yMax, fpc_size)

            fpcGridding(row, col, hag, weight, fpcFP, canopyThreshold)
            canopyWeight = fpcFP[0]
            totalWeight = fpcFP[1]
            nullArray = totalWeight == 0
            totalWeight[nullArray] = 1
            fpc = 100 * (canopyWeight / totalWeight)
            fpc[nullArray] = nullVal

            fpc_flights[f] = fpc

        # Calculate mean FPC of flights, ignoring no-data pixels
        numFpc = np.sum((fpc_flights != nullVal).astype(np.uint8), axis=0)
        nodata = (fpc_flights == nullVal).astype(np.uint8)
        fpc_flights[nodata == 1] = 0
        sumFpc = np.sum(fpc_flights, axis=0)
        numFpc[numFpc == 0] = 1
        fpc = sumFpc / numFpc
        fpc[np.sum(nodata, axis=0) == numFlights] = nullVal

    else:
        x = data["X"]
        y = data["Y"]

        # Ensure no zero values in NUMBER_OF_RETURNS
        if np.min(data["NUMBER_OF_RETURNS"]) < 1:
            data["NUMBER_OF_RETURNS"] += 1

        numberOfReturns = data["NUMBER_OF_RETURNS"]
        weight = 1 / numberOfReturns.astype(np.float32)

        row, col = gridding_methods.xyToRowCol(x, y, xMin, yMax, fpc_size)
        fpc = np.zeros((2, nRows, nCols), dtype=np.float32)
        fpcGridding(row, col, heightAboveGround, weight, fpc, canopyThreshold)

        canopyWeight = fpc[0]
        totalWeight = fpc[1]
        nullArray = totalWeight == 0
        totalWeight[nullArray] = 1
        fpc = 100 * (canopyWeight / totalWeight)
        fpc[nullArray] = nullVal

    return fpc


def makeFPCcolorTable():
    """
    Create an array of the colours we use for the FPC images.

    The returned array has shape (256, 3), and type uint8. The columns are
    red/green/blue values, in the range [0, 255]. This is suitable for use
    with the rw_image_methods.setColorTable function.

    Returns:
      clrTblArr: Array of RGB values
    """
    clrTblArr = np.full((256, 3), 255, dtype=np.uint8)
    clrTblArr[10:90, 0] = np.mgrid[255:0:-80j].round().astype(np.uint8)
    clrTblArr[90:101, 0] = 0
    clrTblArr[:, 2] = clrTblArr[:, 0]
    clrTblArr[10:90, 1] = np.mgrid[255:100:-80j].round().astype(np.uint8)
    clrTblArr[90:101, 1] = 100
    clrTblArr[0, :] = [210, 180, 140]  # Add brown for zero FPC

    return clrTblArr
