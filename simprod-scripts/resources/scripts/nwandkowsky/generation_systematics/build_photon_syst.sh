#!/bin/sh

jobname='NuE'

f='NuE'
e='low_energy'
gcd='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz'
icemodel='spice3_2'
name=''
domeff=0.99
holeice='flasher_p1=0.3_p2=0.0'

if [ $f == 'NuE' ] && [ $e == 'low_energy' ]; then
    seed=22
    domos=1.
    numevents=50000
    gamma=2
    energy_min=0.1
    energy_max=5.
    memory=1000.
fi

if [ $f == 'NuE' ] && [ $e == 'medium_energy' ]; then
    seed=2
    domos=5.
    numevents=10000
    gamma=1.5
    energy_min=5.
    energy_max=10000.
    memory=5000.
fi

if [ $f == 'NuE' ] && [ $e == 'high_energy' ]; then
    seed=222
    domos=5.
    numevents=10
    gamma=1.
    energy_min=1000.
    energy_max=100000.
    memory=6000.
fi

if [ $f == 'NuMu' ] && [ $e == 'low_energy' ]; then
    seed=11
    domos=1.
    numevents=100000
    gamma=2
    energy_min=0.1
    energy_max=5.
    memory=1000.
fi

if [ $f == 'NuMu' ] && [ $e == 'medium_energy' ]; then
    seed=1
    domos=5.
    numevents=20000
    gamma=1.5
    energy_min=5.
    energy_max=10000.
    memory=5000.
fi

if [ $f == 'NuTau' ] && [ $e == 'low_energy' ]; then
    seed=33
    domos=1.
    numevents=50000
    gamma=2
    energy_min=0.1
    energy_max=5.
    memory=1000.
fi

if [ $f == 'NuTau' ] && [ $e == 'medium_energy' ]; then
    seed=3
    domos=5.
    numevents=10000
    gamma=1.5
    energy_min=5.
    energy_max=10000.
    memory=5000.
fi

#k=()

#for i in "${k[@]}"
for i in {7001..8000}
	do
		((j=$i/1000+1))
		JOBID_P=$jobname$i
		echo JOB $JOBID_P photon_syst.condor
		infile="/data/ana/Cscd/StartingEvents/NuGen_new/$f/$e/photon_spice3_2/$j/photon_`printf %08d $i`.i3.zst"
		phofile="/data/ana/Cscd/StartingEvents/NuGen_new/$f/$e/photon_$name/$j/photon_`printf %08d $i`.i3.zst"
		echo VARS $JOBID_P hi=\"$holeice\" memory=\"$memory\" numevents=\"$numevents\" gamma=\"$gamma\" emin=\"$energy_min\" emax=\"$energy_max\" domos=\"$domos\" flavor=\"$f\" seed=\"$seed\" runnumber=\"$i\" outfile=\"$phofile\" domeff=\"$domeff\" ice=\"$icemodel\" gcd=\"$gcd\" infile=\"$infile\" 
                echo Retry $JOBID_P 1
	done
