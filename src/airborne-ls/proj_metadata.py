#!/usr/bin/env python

"""
Not complete

how/where to manage metadata 


"""

import os
import tomllib
import sys
import argparse

from pathlib import Path 
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.ERROR,
                    format='%(asctime)s: %(name)20s: %(levelname)10s: %(message)s')


intColumns = set([x.lower() for x in ['DayStart', 'DayEnd', 'MonthStart', 'MonthEnd', 'YearStart', 'YearEnd',
    'Zone', 'NumRawLasFiles', 'ProvidedEPSG']])
realColumns = set([x.lower() for x in ['OrigTileSize_km', 'FlyingHeight_m', 'SideOverlap_percent', 'FieldOfView_deg',
    'SwathWidth_m', 'AvgPtDensity_pm2', 'LasVersion', 'NonGroundCutoff_m', 'VerticalAccuracy_m',
    'VerticalConfidenceInterval_percent', 'HorizontalAccuracy_m',
    'HorizontalConfidenceInterval_percent', 'FootprintSize_m', 'ScanRate_Hz',
    'PulseRateFrequency_kHz', ]])
textColumns = set([x.lower() for x in ['Project', 'GeneralProjectName_Folder', 'SpecificProjectName',
    'Provider', 'Sensor',
    'SensorCode', 'HorizontalProjection', 'HorizontalDatum', 'VerticalDatum', 'GeoidModel',
    'IMU', 'Classified',  'Thinned',     
    'DataCanBeShared', 'ProjectNotes']])




def build_metadata_template():
   
 
    """
    numerous ways to do this.. good to have input / consensus / discussions 

    write out template for metadata entry - do not use commas within text

    incomplete 

    """
    clist='' # build template 

    outfile=Path(cmdargs.indir).joinpath("%s_metadata.csv" % (cmdargs.proj_name))
    f = open(outfile, 'w')
    for col in clist:
        f.write("%s, None\n" % col[0])
    f.close()


