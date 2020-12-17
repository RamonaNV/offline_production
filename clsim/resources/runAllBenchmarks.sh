#!/bin/sh

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