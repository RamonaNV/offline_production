#!/usr/bin/env python

from icecube import icetray, dataio, dataclasses, phys_services, corsika_reader
from os.path import expandvars


if __name__ == "__main__":
    from optparse import OptionParser
    from I3Tray import *
    from icecube import dataio, corsika_reader,phys_services

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)

    parser.add_option("-o", "--outfile",
                  default="test_flashes.i3", 
                  dest="OUTFILE", 
                  help="Write output to OUTFILE (.i3{.gz} format)")
    parser.add_option("-i", "--infile",
                  default=expandvars('$I3_TESTDATA/DAT010000'),
                  dest="INFILE", 
                  help="Read input from INFILE (.i3{.gz} format)")
    parser.add_option("-r", "--runid",
                  type=int,
                  dest="RUNID", 
                  help="Run Id")
    parser.add_option("-g", "--gcd",
                  default=expandvars("$I3_TESTDATA/sim/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz"),
                  dest="GCDFILE", 
                  help="Read geometry from GCDFILE (.i3{.gz} format)")

    parser.add_option("-s", "--oversample",
                  default=1,
                  type=int,
                  dest="OVERSAMPLE", 
                  help="Oversample showers randomly over the cylinder surface")

    parser.add_option("-n", "--nevents",
                  default=1,
                  type=int,
                  dest="NEVENTS", 
                  help="Number of events to read.")



    (options,args) = parser.parse_args()
    tray = I3Tray()

    tray.context['I3RandomService'] = phys_services.I3GSLRandomService(42)

    tray.Add(corsika_reader.ReadCorsika, "reader", GCDFILE=options.GCDFILE, FileNameList=[options.INFILE], 
             OverSampling=options.OVERSAMPLE, NEvents=options.NEVENTS)

    tray.AddModule("Dump")

    def primary_info(fr):
        print(dataclasses.get_most_energetic_primary(fr['I3MCTree'])) 

    tray.AddModule(primary_info, Streams=[icetray.I3Frame.DAQ])
    tray.AddModule("I3Writer", Filename = options.OUTFILE)
    tray.Execute()
    tray.Finish()
