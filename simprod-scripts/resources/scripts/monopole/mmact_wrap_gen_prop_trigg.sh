#!/bin/bash
# This is a wrapper script for submitting the MMACT generation script (necessary for the OGPU parameter)

while getopts n:d:p:r:g:m:v: option
do
case "${option}"
in
n) numevents=$OPTARG;;
d) outputdir=$OPTARG;;
p) iprocess=$OPTARG;;
r) nprocess=$OPTARG;;
g) gcdfilepath=$OPTARG;;
m) icemodelpath=$OPTARG;;
v) verbosity=$OPTARG;;
esac
done

export OGPU=1

export exestring="mmact_gen_prop_trigg.py"
if [ ${#numevents}    != 0 ] ; then export exestring="$exestring -n $numevents"    ; fi
if [ ${#outputdir}    != 0 ] ; then export exestring="$exestring -d $outputdir"    ; fi
if [ ${#iprocess}     != 0 ] ; then export exestring="$exestring -p $iprocess"     ; fi
if [ ${#nprocess}     != 0 ] ; then export exestring="$exestring -r $nprocess"     ; fi
if [ ${#gcdfilepath}  != 0 ] ; then export exestring="$exestring -g $gcdfilepath"  ; fi
if [ ${#icemodelpath} != 0 ] ; then export exestring="$exestring -m $icemodelpath" ; fi
if [ ${#verbosity}    != 0 ] ; then export exestring="$exestring -v $verbosity"    ; fi


echo -e " ------------\n| MMACT WRAP | About to run the following script\n|            |     $exestring\n ------------\n"
$exestring
export WRAPPED_SCRIPT_EXIT_STATUS=$?
echo -e "\n ------------\n| MMACT WRAP | Finished running the following script with status $WRAPPED_SCRIPT_EXIT_STATUS:\n|            |     $exestring\n ------------"

