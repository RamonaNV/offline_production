#!/bin/sh

# /*The MIT License (MIT)

# Copyright (c) 2020, Hendrik Schwanekamp hschwanekamp@nvidia.com, Ramona Hohl rhohl@nvidia.com

# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGSEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# */

nt=10

while getopts j: flag
do
    case "${flag}" in
        j) nt=${OPTARG};;
        ?) echo "use -j to specify thread count for compiling" && exit 1
    esac
done

echo "Using: $nt threads to compile.";
mkdir -p benchmark


cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=ON -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/opencl_default.txt

cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=ON -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/opencl_shuffled.txt

cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=ON \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/opencl_shuffled32.txt


cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/noJq_noBlocking_default.txt

cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=ON -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/noJq_noBlocking_shuffled.txt

cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=ON ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/noJq_noBlocking_shuffled32.txt


cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=ON \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/noJq_blocking_default.txt

cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=ON \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=ON -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/noJq_blocking_shuffled.txt

cmake -DUSE_JOBQUEUE=OFF -DBLOCK_RANDOM_NUMBERS=ON \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=ON ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/noJq_blocking_shuffled32.txt


cmake -DUSE_JOBQUEUE=ON -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/jq_noBlocking_default.txt

cmake -DUSE_JOBQUEUE=ON -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=ON -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/jq_noBlocking_shuffled.txt

cmake -DUSE_JOBQUEUE=ON -DBLOCK_RANDOM_NUMBERS=OFF \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=ON ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/jq_noBlocking_shuffled32.txt


cmake -DUSE_JOBQUEUE=ON -DBLOCK_RANDOM_NUMBERS=ON \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/jq_blocking_default.txt

cmake -DUSE_JOBQUEUE=ON -DBLOCK_RANDOM_NUMBERS=ON \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=ON -DBENCHMARK_SHUFFLED_32=OFF ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/jq_blocking_shuffled.txt

cmake -DUSE_JOBQUEUE=ON -DBLOCK_RANDOM_NUMBERS=ON \
    -DBENCHMARK_OPENCL=OFF -DBENCHMARK_OPENCL_SHUFFLED=OFF -DBENCHMARK_OPENCL_SHUFFLED_32=OFF \
    -DBENCHMARK_SHUFFLED=OFF -DBENCHMARK_SHUFFLED_32=ON ..
make -j ${nt} 
python clsim/resources/runSmallSim.py --energy=1e4 -r=5 > benchmark/jq_blocking_shuffled32.txt

echo "All benchmarks complete!"