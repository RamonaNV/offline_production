import sys
from I3Tray import I3Units
from icecube import dataclasses as dc

def StressTestHitGenerator( frame, \
                            hit_times = None, weights = None, ) :

    # if weights not set, set it to list filled with None with same size as hit_times
    weights = weights if weights else  [ None for h in hit_times ] 
     
    if len( hit_times ) != len( weights ) :
        print("ERROR : Different sizes for 'hit_times' and 'weights'.")
        print("Either 'weights' needs to be empty (unweighted) or of the same size as 'hit_times'.")
        sys.exit(0)

    # this goes in the frame
    mchit_map = dc.I3MCHitSeriesMap()

    mchits = dc.I3MCHitSeries()
    for t,w in zip( hit_times, weights ) :
        hit = dc.I3MCHit()
        hit.time = t
        hit.npe = w if w else 1
        mchits.append( hit )

    for omkey, geo in frame.Get("I3Geometry").omgeo :
        mchit_map[ omkey ] = mchits

    # calculate bin width
    if weights[0] :
        DT = hit_times[1] - hit_times[0]
        for t0,t1 in zip( hit_times[:-1], hit_times[1:]) :
            if abs( t1 - t0 - DT ) > 1 * I3Units.nanosecond :
                print("ERROR : variable length binning is not supported in simulation")
                sys.exit()

        frame["HitBinWidth"] = dc.I3Double( DT )                

    frame["MCHitSeriesMap"] = mchit_map
    
