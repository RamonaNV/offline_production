In this directory you can find three files (among others):

  -> mmact_gen_prop_trigg.py: This is where the monopoles are generated,
                              propagated and triggered (as the name suggests)
  -> mmact_utils.py: This contains some helper functions and default settings
  -> mmact_wrap_gen_prop_trigg.sh: This is the script you should run - it is a
                                   wrapper around mmact_gen_prop_trigg.py

The purpose of mmact_wrap_gen_prop_trigg.sh is to set the OGPU environment
variable that is required for PPC, and forwards all input arguments to
mmact_gen_prop_trigg.py, which then runs the I3Tray.

The input arguments to mmact_wrap_gen_prop_trigg.sh follow:
  -n: numevents - the number of events that should be run (MUST BE INPUT)
  -d: outputdir - the directory where the output I3 files should be produced
      (MUST BE INPUT)
  -p: iprocess - the number of the current process, required for the RNG (MUST
      BE INPUT)
  -r: nprocess - the total number of processes running, required for the RNG
      (MUST BE INPUT)
  -g: GCD file - the path to the GCD file (default-value: "/data/sim/sim-new/
      downloads/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz")
  -m: ice-model - the path to the ice-model (default-value: "$I3_BUILD/ppc/
      resources/ice/")
  -v: verbosity - currently 0 means print nothing while 1 and 2 mean print all
      (default value: 1)
