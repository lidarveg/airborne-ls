# alsv_addhag
Adds to the given LAZ file(s) to record the estimated height above ground (HAG) for every point.
All points which are classified as "ground" are given HAG of zero. 

Ground points are separated from non-ground points, and used to interpolate an estimate of
ground elevation at the exact (x, y) location of every non-ground point. This is then
subtracted from the elevation of the non-ground point to estimate the height above ground for
each point.

The HAG values for every point are then packed into an EVLR record, and the whole file is
re-written with the points as they were and the extra EVLR record. The values
are scled to whole centimetres and packed as int16 values. The user_id & record_id for
this record are defined in `airborne_ls.const`.

Any points for which a ground elevation value could not be interpolated have their HAG value set
to zero. We should possibly have a missing value indicator instead, but so far have not
done so.

All interpolations are done using natural neighbour interpolation.

