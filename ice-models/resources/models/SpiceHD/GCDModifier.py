import os
import glob
import sys
from I3Tray import *
from icecube import icetray, dataclasses, dataio, simclasses
import argparse

def geometrywriter(frame, geometry="I3Geometry"):
        
        det_geo = frame.Get(geometry)
        frame.Delete(geometry)
      
        dom_map= det_geo.omgeo
        
        print len(dom_map)
        
        for om,pos in dom_map:
	    if verbose: print "in:  ",om,dom_map[om].position
            for line in open("geo-f2k"):
                content=line.split()
                if not(int(content[5])== int(om[0]) and int(content[6])== int(om[1])):
                    continue
                newx=float(content[2])
                newy=float(content[3])
            dom_map[om].position[0]=newx
            dom_map[om].position[1]=newy
            if verbose: print "out: ",om,dom_map[om].position
            if verbose: print ""
            
        det_geo.omgeo=dom_map
            
        frame.Put(geometry,det_geo);

readme="Modify I3 GCD files for use with SpiceHD detector geometry"
parser = argparse.ArgumentParser(description=readme)
parser.add_argument('-i','--inputGCD' ,type=str, required=True,help="name (and location) of input GCD file")
parser.add_argument('-o','--outputGCD',type=str, required=True,help="name (and location) of output GCD file")
parser.add_argument('-v','--verbose',type=bool,default=False)
args = parser.parse_args()

infileGCD=args.inputGCD
outGCD=args.outputGCD
verbose=args.verbose

tray = I3Tray()
tray.AddModule('I3Reader', 'reader',  Filename =infileGCD)
tray.AddModule(geometrywriter,"geometrywriter",Streams=[icetray.I3Frame.Geometry],geometry="I3Geometry")
tray.AddModule("I3Writer", "writer",Streams=[icetray.I3Frame.Geometry, icetray.I3Frame.Calibration,icetray.I3Frame.DetectorStatus],FileName = outGCD)

tray.Execute()

