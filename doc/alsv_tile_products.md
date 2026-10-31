# alsv_tile_products
This script takes the standardised LAZ files (including the estimated height above ground) and
produces a number of raster products.

## DEM
One of the major products is a Digital Elevation Model (DEM). This is interpolated from the
(X, Y, Z) locations of all points classified in the LAZ file as "ground" points.
It is sampled at the requested pixel size, using Natural Neighbour interpolation. This is usually
calculated at fairly high resolution, with the default pixel size being 0.5 metres.

## Maximum Height
This is simply the maximum elevation (Z coordinate) (metres) of any point within the pixel.
Note that these are elevations, not height above the ground.

## Intensity
This is the point return intensity for the point with the maximum height.

## First Return Density
For each pixel, the number of first returns (i.e. points with return number of 1) which
fall within that pixel.

## Ground Returns
This is a boolean raster, showing whether each pixel has one or more ground returns. If
a pixel is False, then no ground returns were within that pixel.

## Non-ground Classes
For each pixel, the most common non-ground class of point returns within that pixel.
A colour table is attached to the image to give the different classes somewhat intuitive colours.

## Canopy Surface model
This is a simple interpolation of the maximum height. The value at each pixel is the natural
neighour interpolated estimate of the maximum height (Z coordinate) for the pixel centre. This
tends to be somewhat smoother and more complete than the maximum height raster described above.

## Foliage Projective Cover
This is calculated using the algorithm of Fisher et al (2020). It uses lidar point returns
to estimate Foliage projective Cover for each pixel, and requires a relatively coarse
pixel size (default is 10 metres). The pixel values are expressed as percentages (i.e.
the percentage of the area which is covered by foliage, projected downwards), and a colour
table is added to identify areas of denser foliage as darker green.

Fisher, A., Armston, J., Goodwin, N., Scarth, P. (2020). *Modelling canopy gap probability,
    foliage projective cover and crown projective cover from airborne lidar metrics in
    Australian forests and woodlands*. Remote Sensing of Environment, 237, 111520.
    doi: 10.1016/j.rse.2019.111520

## Vegetation Height Percentiles
This product is intended to show the distribution of vegetation canopy by height above ground,
at each pixel.

The n-th height percentile is the height above ground which is above <n> percent of the
point returns. This is calculated for a selected set of percentages (1, 5, 25, 50, 75, 95, 99).
The pixel size is generally larger than for products above, with a default of 5 metres.

Point heights below 0.5m are ignored for the purposes of this calculation. Apart from that,
all point returns are included, regardless of their point classification. Pixels with
fewer than 4 returns are assigned the null value.

The output is a raster for each percentile value, with the pixel values being the height
for that percentile, on that pixel.
