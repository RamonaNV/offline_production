# offline_production
Slim version of IceCube's physics production software, including only the projects required to run up to photon generation. 

Full documentation of all projects are available at https://docs.icecube.aq/trunk .

# Dependencies
* cmake 
* python3
* numpy
* boost
* zmq
* gsl
* OpenCL 

## Tested on Ubuntu 20.04

```sh
apt-get install  build-essential cmake libbz2-dev libgsl0-dev libcfitsio-dev
    libboost-system-dev libboost-thread-dev libboost-date-time-dev libzmq5-dev
    libboost-python-dev libboost-serialization-dev libboost-filesystem-dev 
    libboost-program-options-dev libboost-regex-dev libboost-iostreams-dev
    opencl-dev python3-numpy sprng2-dev libarchive-dev
```

## Building

```sh
  mkdir build
  cd build
  export I3_TESTDATA=<path_to_source>/test-data
  cmake <path_to_source>
  make
```

## Running benchmarks

```sh
  ./env-shell.sh
  ./clsim/resources/scripts/benchmark.py --gcd-file=$I3_TESTDATA/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz --numevents 100
```
This should produce output similar to this:

```sh
WARN (clsim): Propagating muons and photons in the same process. This may starve your GPU. (I3CLSimMakePhotons.py:315 in I3CLSimMakePhotons)
 
# these numbers are performance figures for the GPU:
time per photon (GPU): 17.113565755277392 ns
photons per second (GPU): 58433176.01369108 photons per second
 
# these numbers include the host utilization and are probably not meaningful for --numevents=1 (the default). You need more events to even out the startup/setup time.
(avg)   time per photon (actual, including under-utilization): 17.276601912161702 ns
(avg)   photons per second (actual, including under-utilization): 57881752.73611296 photons per second
(total) host time: 12.2 s
(total) waiting time: 91.9 µs (0.001%)
(total) number of kernel calls: 4
        wallclock time: 12.9 s
(avg)   device utilization: 99.05368708262631 %
```
