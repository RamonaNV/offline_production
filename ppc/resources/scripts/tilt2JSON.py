#!/usr/bin/env python

import json
from optparse import OptionParser
import sys
import math


def getKey(tokenList):
    return "%02d,%02d" % (int(tokenList[0]), int(tokenList[1]))


def main(outFile, tiltFile, cableFile):

    tiltMap = {}
    with open(tiltFile) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) != 6:
                print "Bad line in %s: %s" % (tiltFile, line)
                sys.exit(-1)
            key = getKey(tokens)
            tiltMap[key] = {}
            # IN GCD database, PMT is at top, so reverse orientation
            tiltMap[key]["nx"] = -1. * float(tokens[2])
            tiltMap[key]["ny"] = -1. * float(tokens[3])
            tiltMap[key]["nz"] = -1. * float(tokens[4])
            tiltMap[key]["tiltUncertainty"] = float(tokens[5]) * math.pi / 180.

    with open(cableFile) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) != 4:
                print "Bad line in %s: %s" % (cableFile, line)
                sys.exit(-1)
            key = getKey(tokens)
            if key not in tiltMap:
                tiltMap[key] = {}
            # DOM x-axis is LED7.  Current analysis reports the
            # cable azimuth, which is assumed to be +90 degrees
            # from LED7.  Input is in degrees.
            azimuth = (float(tokens[2]) - 90) * math.pi / 180.
            # Report -pi < azimuth < pi
            azimuth = math.atan2(math.sin(azimuth), math.cos(azimuth))
            azimuthUncertainty = float(tokens[3]) * math.pi / 180.
            tiltMap[key]["azimuth"] = azimuth
            tiltMap[key]["azimuthUncertainty"] = azimuthUncertainty

    with open(outFile, "w") as out:
        out.write(json.dumps(tiltMap, sort_keys=True,
                  indent=4, separators=(',', ': ')))


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-o", dest="output", help="Output JSON file name")
    parser.add_option("-t", dest="tilt", help="DOM tilt file")
    parser.add_option("-c", dest="cable", help="Cable azimuth file")
    (options, args) = parser.parse_args()
    if any(x is None for x in [options.output, options.tilt, options.cable]):
        parser.print_help()
        sys.exit(-1)
    main(options.output, options.tilt, options.cable)
