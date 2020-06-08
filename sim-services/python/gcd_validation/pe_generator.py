import sys
from I3Tray import I3Units
from icecube import simclasses as sc

def StressTestPEGenerator( frame, \
                            hit_times = None, weights = None, ) :

    # if weights not set, set it to list filled with None with same size as hit_times
    weights = weights if weights else  [ None for h in hit_times ] 
     
    if len( hit_times ) != len( weights ) :
        print("ERROR : Different sizes for 'hit_times' and 'weights'.")
        print("Either 'weights' needs to be empty (unweighted) or of the same size as 'hit_times'.")
        sys.exit(0)

    mcpes = sc.I3MCPESeries()
    for t,w in zip( hit_times, weights ) :
        hit = sc.I3MCPE()
        hit.time = t
        hit.npe = w if w else 1
        mcpes.append( hit )

    # this goes in the frame
    mcpe_map = sc.I3MCPESeriesMap()
    for omkey, geo in frame.Get("I3Geometry").omgeo :
        # Only InIce and IceTop DOMs
        if 0 < omkey.om < 65:
            mcpe_map[ omkey ] = mcpes

    frame["I3MCPESeriesMap"] = mcpe_map
    
