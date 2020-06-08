#!/bin/sh

jobname_D='D'

f='NuE'
s=22
e='low_energy'
y='IC86_flasher_p1=0.3_p2=1.0'
p='photon_spice3_2'
phot='mcpe'
domeff=0.99
domos=1
holeice='flasher_p1=0.3_p2=1.0'

for i in {1..8000}
	do
	    ((j=$i/1000+1))
	    ((l=$((i-$((j-1))*1000))))
	    JOBID_D=$jobname_D$i
            echo JOB $JOBID_D detector_mcpe.condor
            phofile="/data/ana/Cscd/StartingEvents/NuGen_new/$f/$e/$p/$j/photon_`printf %08d $i`.i3.zst"
            detfile="/data/ana/Cscd/StartingEvents/NuGen_new/$f/$e/$y/detector/$j/det_`printf %08d $i`.i3.zst"
	    bgfile="/data/ana/Cscd/StartingEvents/CORSIKA_bg/12531/"$phot"_spice3_2/1/"$phot"_0000`printf %04d $l`.i3.zst"
            echo VARS $JOBID_D seed=\"$s\" runnumber=\"$i\" infile=\"$phofile\" outfile=\"$detfile\" bgfile=\"$bgfile\" hi=\"$holeice\" domeff=\"$domeff\" domos=\"$domos\"
            echo Retry $JOBID_D 2

        done
