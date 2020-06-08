#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Updates official GCD checksum file.')
parser.add_argument('gcd_filename', help='Name of GCD file.')
args = parser.parse_args()

import os
import urllib
import pickle
from icecube import icetray
from icecube import dataio

# get the checksum file
url = 'http://simprod.icecube.wisc.edu/downloads/GCD/checksums'
response = urllib.urlopen(url).read()
title = str(response).split('<title>')[1].split('</title>')[0]
checksums = pickle.loads(response) if '404' not in title\
            else {}

gcd = dataio.I3File(args.gcd_filename)

geometry_frame = None
calibration_frame = None
detector_status_frame = None

while gcd.more():
    frame = gcd.pop_frame()
    if 'I3Geometry' in frame:
        geometry_frame = frame['I3Geometry']
    if 'I3Calibration' in frame:
        calibration_frame = frame['I3Calibration']
    if 'I3DetectorStatus' in frame:
        detector_status_frame = frame['I3DetectorStatus']

#m = hashlib.sha256()
geometry_hash = geometry_frame.sha256()
calibration_hash = calibration_frame.sha256()
detector_status_hash = detector_status_frame.sha256()

print("geo hash = %s" % geometry_hash )
print("cal hash = %s" % calibration_hash )
print("detstat hash = %s" % detector_status_hash )



