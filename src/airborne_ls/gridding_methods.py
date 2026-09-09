#!/usr/bin/env python

"""
This module stores a series of functions for converting point clouds into gridded products.
"""

import numpy as np
import pynninterp
from numba import njit
from scipy import ndimage

from airborne_ls import const

# ----------------------------------------------------------------------------------------------------
# DEM Generation Functions
# ----------------------------------------------------------------------------------------------------


def makeDem(xVals, yVals, zVals, header, tile_s, psize, nullVal=-999.0):
    """
    Interpolate the ground returns on a regular grid to make a DEM image.
    """
    nRows = int(np.ceil(tile_s / psize))
    pxlCoords = get_grid(nRows, header.x_min, header.y_max + 1, psize)
    dem = pynninterp.NaturalNeighbour(xVals, yVals, zVals,
        pxlCoords[0].astype(np.float64), pxlCoords[1].astype(np.float64))
    dem[np.isnan(dem)] = nullVal
    return dem


def makeDemTile(xVals, yVals, zVals, x_min, y_max, tile_s, psize,
        nullVal=-999.0, Linear=False):
    """
    Interpolate the ground returns on a regular grid to make a DEM image.
    """
    nRows = int(np.ceil(tile_s / psize))
    pxlCoords = get_grid(nRows, x_min, y_max, psize)

    # Rounding to fix NaN issues
    xVals = np.round(xVals, 2)
    yVals = np.round(yVals, 2)
    zVals = np.round(zVals, 2)

    if Linear:
        dem = pynninterp.Linear(xVals, yVals, zVals,
            pxlCoords[0].astype(np.float64), pxlCoords[1].astype(np.float64))
    else:
        dem = pynninterp.NaturalNeighbour(xVals, yVals, zVals,
            pxlCoords[0].astype(np.float64), pxlCoords[1].astype(np.float64))

    dem[np.isnan(dem)] = nullVal
    return dem


def interpPoints(xVals, yVals, zVals, x_min, y_max, tile_s, psize, nullVal=-999.0):
    """
    Interpolate the ground returns on a regular grid to make a DEM image.
    """
    nRows = int(np.ceil(tile_s / psize))
    pxlCoords = get_grid(nRows, x_min, y_max + 1, psize)

    dem = pynninterp.NaturalNeighbour(xVals, yVals, zVals,
        pxlCoords[0].astype(np.float64), pxlCoords[1].astype(np.float64))

    dem[np.isnan(dem)] = nullVal
    return dem


# ----------------------------------------------------------------------------------------------------
# Grid and Coordinate Conversion Functions
# ----------------------------------------------------------------------------------------------------


def xyToRowCol(x, y, xMin, yMax, pixSize):
    """
    Convert arrays of x and y into arrays of row and column in a regular grid.
    """
    col = (np.floor((x - xMin) / pixSize)).astype(np.uint32)
    row = (np.floor((yMax - y) / pixSize)).astype(np.uint32)
    return row, col


def get_grid(xr, xst, yst, psize):
    """
    Build paired x and y grid locations for interpolation with pynninterp.
    """
    yr = xr
    x_id = np.arange(int(xr)) * psize + xst
    x_ids = (0.5 * psize) + np.ones((int(xr), 1), np.float64) * x_id
    y_id = np.floor(yst) - np.arange(int(yr)) * psize
    y_ids = np.ones((int(yr), 1), np.float64) * y_id - (0.5 * psize)
    return [x_ids, np.rot90(y_ids, 3)]


# ----------------------------------------------------------------------------------------------------
# Filtering Functions
# ----------------------------------------------------------------------------------------------------


@njit
def runPitInfill(ny, nx, ksize, med):
    """
    Fill pits in the DEM using a neighbourhood-based approach.
    """
    for p in range(ksize, ny - ksize, ksize):
        for q in range(ksize, nx - ksize, ksize):
            k = med[p - ksize : p + (ksize + 1), q - ksize : q + (ksize + 1)]
            if np.sum(np.isnan(k)) < 6 and k[ksize, ksize] == np.nanmin(k):
                k[ksize, ksize] = 99999.0
                med[p, q] = np.nanmin(k) + 0.01


def runDEM_filter(ny, nx, ksize, dem):
    """
    Smooth the DEM to reduce roughness for better flow direction.
    """
    avDem = np.zeros((ny + 2 * ksize, nx + 2 * ksize), dtype=np.float32)
    ct = np.zeros((ny + 2 * ksize, nx + 2 * ksize), dtype=np.float32)
    dt = (dem > 0).astype(np.float32)
    kernel = circleLocs(ksize)

    for p in range(2 * ksize + 1):
        for q in range(2 * ksize + 1):
            if kernel[p, q]:
                avDem[p : ny + p, q : nx + q] += dem
                ct[p : ny + p, q : nx + q] += dt

    avDem = avDem / ct
    avDem = avDem[ksize : (ny + ksize), ksize : (nx + ksize)]
    return avDem


def circleLocs(radius):
    """
    Generate a circular kernel for neighbourhood operations.
    """
    y, x = np.ogrid[-radius : radius + 1, -radius : radius + 1]
    return x**2 + y**2 <= radius**2


###################################################################################################
@njit
def maxH_workflow_layers(row, col, x, y, z, intensity, classi,
        xArr, yArr, zArr, intensityAtMaxH, haveGroundReturn):
    """
    Compute the maximum height grid and associated x, y locations from a LiDAR point cloud as
    2D arrays.

    Parameters:
        row, col: Arrays of row and column indices for each point.
        x, y, z: Arrays of x, y, and z coordinates of the LiDAR points.
        intensity: Array of intensity values for each point.
        classi: Array of classification values for each point.
        xArr, yArr, zArr: 2D arrays to store x, y, and z values at the maximum height for each
                          grid cell.
        intensityAtMaxH: 2D array to store intensity values at the maximum height for each
                         grid cell.
        haveGroundReturn: 2D array to indicate whether a ground return exists for each grid cell.

    Returns:
        Updates the input arrays in place with the computed values.

    """
    numPts = len(row)
    for i in range(numPts):
        (r, c) = (row[i], col[i])
        if z[i] > zArr[r, c]:
            xArr[r, c] = x[i]
            yArr[r, c] = y[i]
            zArr[r, c] = z[i]
            intensityAtMaxH[r, c] = intensity[i]
        if classi[i] == const.PTCLASS_GROUND:
            haveGroundReturn[r, c] = 1


@njit
def makeClassCounts(row, col, x, y, z, pntClass, classCounts):
    """
    Make a per-pixel counts of point classification values

    Assumes that classCounts has been intialized to all zeros, and updates
    it in-place.

    Parameters:
      row, col: Arrays of the pixel row and column values for each point
                return (numPoints)
      x, y, z: Arrays of the coordinates of each point return (numPoints)
      pntClass: Array of classification value for each point return (numPoints)
      classCounts: Array of per-pixel counts for each classification value,
                   shape (nClass, nRows, nCols).
    """
    classUpperBound = classCounts.shape[0]

    numPts = len(row)
    for i in range(numPts):
        (r, c) = (row[i], col[i])
        cls = pntClass[i]
        if cls < classUpperBound:
            classCounts[cls, r, c] += 1


@njit
def fstR_density(row, col, x, y, density):
    """
    Increment the density count for each (row, col) grid cell based on the first return data.

    Parameters:
        row, col: Arrays of row and column indices for each point.
        x, y: Arrays of x and y coordinates of the first return data points
              (not used in this function).
        density: 2D array to store the density count for each grid cell.
    """
    numPts = len(row)
    for i in range(numPts):
        r, c = row[i], col[i]
        density[r, c] += 1


###################################################################################################
def createHeightAboveGround(nonGround, x, y, z, xVals, yVals, zVals):
    """
    Create height above ground by interpolating ground elevation for the coordinates
    of each non-ground point and calculating the difference.

    Parameters:
        nonGround: Boolean array indicating non-ground points.
        x, y, z: Arrays of x, y, and z coordinates of all points.
        xVals, yVals, zVals: Arrays of x, y, and z coordinates of ground points.

    Returns:
        heightAboveGround: Array of height above ground for all points.
    """
    # Count the number of non-ground points
    nNonGrd = np.sum(nonGround)

    # Prepare an array of x, y coordinates for non-ground points
    xys = np.zeros((nNonGrd, 2), dtype=np.float64)
    xys[:, 0] = x[nonGround]
    xys[:, 1] = y[nonGround]

    # Extract z values for non-ground points
    ztemp = z[nonGround]

    # Initialise the heightAboveGround array
    heightAboveGround = np.zeros(z.shape, dtype=z.dtype)

    # Process in chunks if the number of non-ground points exceeds the chunk size
    chunkSize = int(5e6)
    if nNonGrd > chunkSize:
        # Calculate the number of chunks and their indices
        nChunks = (nNonGrd // chunkSize) + 1
        incs = np.concatenate([np.arange(0, nChunks) * chunkSize, [nNonGrd]]).astype(
            np.int32
        )

        # Process each chunk
        for i in range(len(incs) - 1):
            start, end = incs[i], incs[i + 1]
            irregZ = pynninterp.NaturalNeighbourPts(
                xVals, yVals, zVals, xys[start:end, :]
            )
            hh = ztemp[start:end] - irregZ
            hh[np.isnan(irregZ)] = 0  # Handle NaN values in interpolated ground heights

            # Concatenate results
            if i == 0:
                heights = hh
            else:
                heights = np.concatenate([heights, hh])

        # Update the heightAboveGround array for non-ground points
        heightAboveGround[nonGround] = heights
    else:
        # Process all points in one go if the number of non-ground points is small
        irregZ = pynninterp.NaturalNeighbourPts(xVals, yVals, zVals, xys)
        ztemp[np.isnan(irregZ)] = 0.0  # Handle NaN values in non-ground heights
        irregZ[np.isnan(irregZ)] = (
            0.0  # Handle NaN values in interpolated ground heights
        )
        heightAboveGround[nonGround] = ztemp - irregZ

    return heightAboveGround


@njit
def maxH_xyzLocs(x, y, z, psize, nullV=-999.0):
    """
    Compute the maximum height (z) and associated x, y locations for each grid cell.

    Parameters:
        x, y, z: Arrays of x, y, and z coordinates of the points.
        psize: Grid cell size.
        nullV: Default value for cells with no data (default is -999.0).

    Returns:
        Tuple of three 1D arrays:
            - xv: x-coordinates of the maximum height points.
            - yv: y-coordinates of the maximum height points.
            - maxH: Maximum height values.
    """
    # Convert x, y coordinates to row and column indices
    xMin = np.floor(x.min())
    yMax = np.ceil(y.max())
    col = (np.floor((x - xMin) / psize)).astype(np.uint32)
    row = (np.floor((yMax - y) / psize)).astype(np.uint32)
    numPts = len(row)

    # Determine the size of the grid
    maxRow, maxCol = np.max(row) + 1, np.max(col) + 1

    # Initialise arrays with null values
    maxH = np.full((maxRow, maxCol), nullV, dtype=np.float32)
    xv = np.full((maxRow, maxCol), nullV, dtype=np.float32)
    yv = np.full((maxRow, maxCol), nullV, dtype=np.float64)

    # Update arrays with maximum height and associated x, y values
    for i in range(numPts):
        r, c = row[i], col[i]
        if z[i] > maxH[r, c]:
            maxH[r, c] = z[i]
            xv[r, c] = x[i]
            yv[r, c] = y[i]

    # Flatten arrays and filter out null values
    xv = xv.flatten()
    yv = yv.flatten()
    maxH = maxH.flatten()
    valid = maxH > nullV

    return (
        xv[valid].astype(np.float64),
        yv[valid].astype(np.float64),
        maxH[valid].astype(np.float64),
    )


@njit
def maxH_array(row, col, z, maxH_hag):
    """
    Update a 2D grid with the maximum height values from a LiDAR point cloud.

    Parameters:
        row: Array of row indices for each point.
        col: Array of column indices for each point.
        z: Array of height (z) values for each point.
        maxH_hag: 2D array to store the maximum height values for each grid cell.

    Updates:
        Modifies the `maxH_hag` array in place by updating each cell with the maximum height value.
    """
    numPts = len(row)
    for i in range(numPts):
        r, c = row[i], col[i]
        maxH_hag[r, c] = max(maxH_hag[r, c], z[i])


###################################################################################################
@njit
def count_fstR(row, col, density):
    """
    Increment the density count for each (row, col) grid cell.

    Parameters:
        row: Array of row indices for each point.
        col: Array of column indices for each point.
        density: 2D array to store the density count for each grid cell.

    Updates:
        Modifies the `density` array in place by incrementing the count for each (row, col) pair.
    """
    for p in range(len(row)):
        density[row[p], col[p]] += 1


###################################################################################################
def doHeightPercentileOutputs(x, y, xMin, yMax, heightAboveGround, tile_s, psize,
        pptiles=(1, 5, 25, 50, 75, 95, 99)):
    """
    Generate gridded outputs of height percentiles.

    Parameters:
        x, y: Arrays of x and y coordinates of the points.
        xMin, yMax: Minimum x and maximum y coordinates of the grid.
        heightAboveGround: Array of height above ground values.
        tile_s: Size of the tile (assumed to be square).
        psize: Grid cell size.
        pptiles: List of percentiles to calculate (default: (1, 5, 25, 50, 75, 95, 99)).

    Returns:
        percentile_arr: 3D array of height percentiles for each grid cell.
    """
    # Calculate the number of rows and columns for the grid
    nRows = int(np.ceil(tile_s / psize))
    nCols = nRows  # Assumes square tiles

    # Convert x, y coordinates to row and column indices
    row, col = xyToRowCol(x, y, xMin, yMax, psize)

    # Initialise the percentile array with a null value
    nullVal = -999
    percentile_arr = np.full((len(pptiles), nRows, nCols), nullVal, dtype=np.float32)

    # Calculate percentiles for each grid cell
    for r in range(nRows):
        thisRow = row == r
        colsThisRow = col[thisRow]
        hgtThisRow = heightAboveGround[thisRow]

        for c in range(nCols):
            thisCol = colsThisRow == c
            hgtThisCell = hgtThisRow[thisCol]
            hgtThisCell = hgtThisCell[hgtThisCell > 0.5]  # Filter heights > 0.5

            if (
                len(hgtThisCell) > 3
            ):  # Only calculate percentiles if there are enough points
                for idx, pp in enumerate(pptiles):
                    percentile_arr[idx, r, c] = np.percentile(hgtThisCell, pp)

    return percentile_arr


@njit
def doHeightPercentileOutputs_idx(x, y, xst_bin, yst_bin, binSize, heightAboveGround,
        ptile_s, pptiles=(1, 5, 25, 50, 75, 95, 99), nullVal=-999.0):
    """
    Generate gridded outputs of height percentiles using bin indices.

    Parameters:
        x, y: Arrays of x and y coordinates of the points.
        xst_bin, yst_bin: Starting x and y coordinates of the bin.
        binSize: Size of the bin (assumed to be square).
        heightAboveGround: Array of height above ground values.
        ptile_s: Grid cell size for percentiles.
        pptiles: List of percentiles to calculate (default: (1, 5, 25, 50, 75, 95, 99)).
        nullVal: Default value for grid cells with no data (default: -999.0).

    Returns:
        percentile_arr: 3D array of height percentiles for each grid cell.
        nRows_pct: Number of rows in the percentile grid.
    """
    # Calculate the number of rows and columns for the percentile grid
    nRows_pct = nCols_pct = int(np.ceil(binSize / ptile_s))

    # Initialise the percentile array with a null value
    percentile_arr = np.full(
        (len(pptiles), nRows_pct, nCols_pct), nullVal, dtype=np.float32
    )

    # Convert x, y coordinates to row and column indices
    col = (np.floor((x - xst_bin) / ptile_s)).astype(np.uint32)
    row = (np.floor((yst_bin - y) / ptile_s)).astype(np.uint32)

    # Reset the percentile array for the given rows and columns
    percentile_arr[:, row, col] = 0

    # Calculate percentiles for each grid cell
    for r in range(nRows_pct):
        thisRow = row == r
        colsThisRow = col[thisRow]
        hgtThisRow = heightAboveGround[thisRow]

        for c in range(nCols_pct):
            thisCol = colsThisRow == c
            hgtThisCell = hgtThisRow[thisCol]
            hgtThisCell = hgtThisCell[hgtThisCell > 0.5]  # Filter heights > 0.5

            if len(hgtThisCell) > 3:
                # Only calculate percentiles if there are enough points
                for idx, pp in enumerate(pptiles):
                    percentile_arr[idx, r, c] = np.percentile(hgtThisCell, pp)

    return percentile_arr, nRows_pct


###################################################################################################
def dem_infill(dem, codes, nullVal=-999.0, minElev=-4):
    """
    Infill missing areas of a Digital Elevation Model (DEM). This is mainly used
    to handle interpolation chunks that are smaller than empty areas (e.g., water bodies).
    Note: It may be more ideal to mask out water bodies manually.

    Parameters:
        dem: 2D or 3D array representing the DEM.
        codes: 2D or 3D array of codes indicating specific areas (e.g., water, noise).
        nullVal: Value to assign to invalid or missing data (default: -999.0).
        minElev: Minimum elevation threshold (default: -4).

    Returns:
        dem: The infilled DEM with missing areas filled.
    """
    # Determine the shape of the DEM
    if len(dem.shape) == 2:
        ny, nx = dem.shape
        nz = 0
    elif len(dem.shape) == 3:
        nz, ny, nx = dem.shape
        dem = dem.reshape(ny, nx)
        codes = codes.reshape(ny, nx)
    else:
        raise ValueError("DEM must be a 2D or 3D array.")

    # Reset problematic code values (e.g., water, noise)
    for code in [7, 9, 18]:
        vals = codes == code
        if np.sum(vals) > 0:
            dem[vals] = nullVal

    # Reset DEM values of 0 to the null value
    vals = dem == 0.0
    if np.sum(vals) > 0:
        dem[vals] = nullVal

    # Identify bad values below the minimum elevation
    bad_vals = dem < minElev
    # Only proceed if there are enough bad values
    if (np.sum(bad_vals) > 3 and np.max(dem) > minElev):
        # Ensure there are valid values to interpolate from
        # Get valid elevation points above the minimum elevation
        valid_indices = np.argwhere(dem > minElev)
        yVals = valid_indices[:, 0].astype(np.float64)
        xVals = valid_indices[:, 1].astype(np.float64)
        zVals = dem[dem > minElev].flatten().astype(np.float64)

        # Generate grid coordinates for interpolation
        nRows = int(np.max([ny, nx]))
        pxlCoords = get_grid(nRows, 0, ny, 1.0)  # Generate linear grid

        # Perform natural neighbour interpolation
        dem2 = pynninterp.Linear(xVals, yVals, zVals,
            pxlCoords[0].astype(np.float64), pxlCoords[1].astype(np.float64))
        dem2 = dem2[0:ny, 0:nx]
        dem2 = np.copy(dem2[::-1, :])  # Flip vertically

        # Update bad values in the DEM with interpolated values
        dem[bad_vals] = dem2[bad_vals]

        # Reset any areas with code 255 to the null value
        vals = codes == 255
        if np.sum(vals) > 0:
            dem[vals] = nullVal

    # Reshape DEM back to its original shape if it was 3D
    if nz > 0:
        dem = dem.reshape(nz, ny, nx)

    return dem


###################################################################################################
# CSM
def chm_alg(chunk, chunk_hag, maxH_hag, psize, xst_bin, yst_bin, binSize, nRows, nCols,
        nullVal, minH_thres=1.0):
    """
    Modified version of the pit-free algorithm for generating a Canopy Height Model (CHM).

    Parameters:
        chunk: Input data chunk containing point cloud data.
        chunk_hag: Height above ground (HAG) values for the chunk.
        maxH_hag: Maximum height above ground for each bin in the chunk.
        psize: Grid cell size.
        xst_bin, yst_bin: Starting x and y coordinates of the bin.
        binSize: Size of the bin (assumed to be square).
        nRows, nCols: Number of rows and columns in the output grid.
        nullVal: Value to assign to empty grid cells.
        minH_thres: Minimum height threshold for processing (default: 1.0).

    Returns:
        outarr: 2D array representing the generated CHM.
    """
    # Assign height above ground values to the 'Z' field in the chunk
    chunk["Z"] = chunk_hag

    # Initialise the output array
    outarr = np.zeros((nRows, nCols), dtype=np.float32)

    # Define the dilation structure based on the grid cell size
    if psize < 0.2:
        dil_struct = np.ones((5, 5))
    else:
        dil_struct = np.ones((4, 4))

    # Determine height increments for processing
    max_height = np.max(chunk["Z"])
    if max_height > 5.0:
        h_incs = np.arange(0, int(np.ceil(max_height / 5) * 5), 5)
    elif max_height > 1.0:
        h_incs = [0]
    else:
        return outarr  # Return empty array if max height is too low

    # Set the first height increment to the minimum height threshold
    h_incs[0] = minH_thres
    nElems = len(chunk)

    # Process each height increment
    for idx, inc in enumerate(h_incs):
        vals = chunk["Z"] >= inc
        if np.sum(vals) > 0 and np.sum(vals) / nElems > 0.05:
            data_sub = chunk[vals]

            # Create a DEM tile for the current height increment
            csm = makeDemTile(data_sub["X"], data_sub["Y"], data_sub["Z"],
                xst_bin, yst_bin, binSize, psize, nullVal=0, Linear=True)

            # Create a mask for areas above the current height increment
            msk = maxH_hag >= inc
            msk = ndimage.binary_dilation(msk * 1, structure=dil_struct)

            # Apply the mask to the current surface model
            csm = csm * msk

            # Update the output array with the maximum values
            outarr[csm > outarr] = csm[csm > outarr]

    return outarr


###################################################################################################
