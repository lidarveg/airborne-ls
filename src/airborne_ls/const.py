"""
Define some constants to be generally available through the package
"""

# Values of the point classification, following the ASPRS standard. See
# https://community.asprs.org/leadership-restricted/leadership-content/public-documents/standards
#
PTCLASS_NEVERCLASSIFIED = 0
PTCLASS_UNCLASSIFIED = 1
PTCLASS_GROUND = 2
PTCLASS_LOWVEGETATION = 3
PTCLASS_MEDIUMVEGETATION = 4
PTCLASS_HIGHVEGETATION = 5
PTCLASS_BUILDING = 6
PTCLASS_NOISE_LOWPOINT = 7
PTCLASS_WATER = 9
PTCLASS_RAIL = 10
PTCLASS_ROADSURFACE = 11
PTCLASS_WIREGUARD = 13
PTCLASS_WIRECONDUCTOR = 14
PTCLASS_TRANSMISSIONTOWER = 15
PTCLASS_WIRESTRUCTURECONNECTOR = 16
PTCLASS_BRIDGEDECK = 17
PTCLASS_NOISE_HIGHPOINT = 18
PTCLASS_OVERHEADSTRUCTURE = 19
PTCLASS_IGNOREDGROUND = 20
PTCLASS_SNOW = 21
PTCLASS_TEMPORALEXCLUSION = 22

# User-defined (i.e. non-standard) classification values
PTCLASS_NOISE_PROVIDERDEFINED = 64


# We define some VLR types, under our own user_id
VLR_USERID_JRSRP = "JRSRP Aust"
# Records for our square bin index.
VLR_RECORDID_BINSIZE = 1            # (float64, uint64)
VLR_RECORDID_BINBOUNDS = 2          # All uint64
