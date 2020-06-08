#!/bin/zsh
#
#$ -S /bin/zsh
#$ -m a
#$ -l arch=amd64
#$ -j y
#$ -R y
#$ -pe multicore 2
# h_vmem is per-core, so divide
#$ -l h_vmem=14G
#$ -l h_stack=20M
#$ -l h_cpu=01:00:00
# -o /afs/ifh.de/group/amanda/scratch/nwhiteho/fittables/log

TABLEDIR=/afs/ifh.de/user/m/middell/scratch/data/photonicstables/AHAv1_level1stub/tables/AHA07v1ice/showers
OUTDIR=/afs/ifh.de/group/amanda/scratch/nwhiteho/fittables-moret-bugfix

#TABLEDIR=/afs/ifh.de/user/m/middell/scratch/data/photonicstables/mdagost_flasher/tables
#TABLE=DOM_37-horiz_only_beamedLED_10-7-405nm_AHA_NOHOLEICE_emit_10cm_receive.pt.norm.prob

cd $TABLEDIR
TABLES=(`ls -1 *.prob | sort`)
#TABLES=( $(for i in `ls -1 $TABLEDIR/*.prob`; do echo $i | xargs basename; done) )

TABLE=${TABLES[$SGE_TASK_ID]}
SPQR_THREADS=1
BLAS_THREADS=$NSLOTS

OUT=$TABLE.pspl.fits

# User servicable options end

export PATH=/afs/ifh.de/group/amanda/software/$OS_ARCH/bin:$PATH
export LD_LIBRARY_PATH=/afs/ifh.de/group/amanda/software/${OS_ARCH}/lib:/afs/ifh.de/user/n/nwhiteho/lib
export PYTHONPATH=/afs/ifh.de/user/n/nwhiteho/python
export GOTO_NUM_THREADS=$BLAS_THREADS
export SPQR_NTHREADS=$SPQR_THREADS

cd /afs/ifh.de/user/n/nwhiteho/photospline/fitter
afscp $TABLEDIR/$TABLE $TMPDIR/input.pt.prob
echo "Fitting table $TABLEDIR/$TABLE..."
(python glam-photonics.py $TMPDIR/input.pt.prob $TMPDIR/output.fits > $TMPDIR/log) || exit 1
afscp $TMPDIR/output.fits $OUTDIR/$OUT && echo "Output moved to $OUTDIR/$OUT."
afscp $TMPDIR/log $OUTDIR/$TABLE.log
