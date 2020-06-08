# offline_production
Slim version of IceCube's physics production software, including only the projects required to run up to photon generation. 

Full documentation of all projects are available at https://docs.icecube.aq .

# Dependencies
* cmake 
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
    opencl-dev
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
  ./clsim/resources/scripts/benchmark.py --gcd-file=$I3_TESTDATA/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz
```
