#!/bin/sh

for s in + -; do
for a in + -; do
dir=err_s${s}.05_a${a}.05
echo $dir
mkdir $dir
cp cfg.txt icemodel.par tilt.dat tilt.par $dir
awk 'BEGIN {o=0; n=0; s=exp('$s'log(1.05)); a=exp('$a'log(1.05)); print s, a}'
awk 'BEGIN {o=0; n=0; s=exp('$s'log(1.05)); a=exp('$a'log(1.05))} {print $1, s*$2, a*($3+o)-o, $4}' icemodel.dat > $dir/icemodel.dat
done
done
