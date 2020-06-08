#!/bin/sh

jobname_L2='L2'

f='NuE'
e='low_energy'
y='IC86_flasher_p1=0.3_p2=1.0'

for i in {1..8000}
	do
	    ((j=$i/1000+1))
            # read the detector files and perform pole processing (filtering) and offline processing (L2 recos)
            JOBID_L2=$jobname_L2$i
            echo JOB $JOBID_L2 l2.condor
	    detfile="/data/ana/Cscd/StartingEvents/NuGen_new/$f/$e/$y/detector/$j/det_`printf %08d $i`.i3.zst"
            l2file="/data/ana/Cscd/StartingEvents/NuGen_new/$f/$e/$y/l2/$j/l2_`printf %08d $i`.i3.zst"
            echo VARS $JOBID_L2 infile=\"$detfile\" outfile=\"$l2file\"  
            echo Retry $JOBID_L2 20

        done
