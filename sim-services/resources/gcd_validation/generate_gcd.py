#!/usr/bin/env python

import sys
try:
    from icecube import I3Db
except ImportError:
    print("ERROR : This script needs I3Db to run.")
    sys.exit(1)

from I3Tray import I3Tray
from I3Tray import I3Units

from icecube import dataclasses as dc

import os
from os.path import expandvars

################################################################################
# 'version' is the latest published version of the GCD.                        #
# So -1 versions mean one hasn't been generated yet.                           #
# https://wiki.icecube.wisc.edu/index.php/GCD_File_Info_for_Production         #
# Only increment the version when you plan to publish a new file.  This way    #
# people can easily recreate the current published file.
################################################################################
season_params = { "2015" : {"MJD" : 57161, "version" : -1 }, \
                  "2014" : {"MJD" : 56784, "version" : -1 }, \
                  "2013" : {"MJD" : 56429, "version" : 1  }, \
                  "2012" : {"MJD" : 56063, "version" : 1  }, \
                  "2011" : {"MJD" : 55697, "version" : 2  }, \
                  "IC79" : {"MJD" : 55380, "version" : -1 }, \
                  "IC59" : {"MJD" : 55000, "version" : -1 }, \
                  "IC40" : {"MJD" : 54649, "version" : -1 } \
                  }

################################################################################
# Get the user input.  The season needs to be one of the keys in the
# season_params dictionary.
################################################################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s","--season",
                  dest="SEASON", default="2014" ,
                  help="Season to generate %s" % season_params.keys())
(options, args) = parser.parse_args()


# Validate the options
if options.SEASON not in season_params :
    msg = "ERROR : Season %s is not supported.  Please choose from the following : \n%s" \
    % (options.SEASON, season_params.keys())
    print(msg)
    sys.exit(1)
 
print("Please be patient.  This will take several minutes...")

MJD = season_params[options.SEASON]["MJD"]
new_version = season_params[options.SEASON]["version"] + 1
out_filename = "GeoCalibDetectorStatus_"+options.SEASON+"."+str(MJD)+"_V"+str(new_version)+".i3.gz"
logfilename = out_filename.replace(".i3.gz",".log")
logfile = open(logfilename, "w")

from icecube import icetray
icetray.logging.rotating_files(logfilename)

import subprocess
SCRIPT_PATH = expandvars("$I3_BUILD/sim-services/resources/gcd_validation/details/")
    
################################################################################
# First get the G, C, and D frames from the DB and generate a raw GCD file.    #
# This also includes the Bad DOM list.                                         #
################################################################################
print("Pulling GCD information from the database...")
script = SCRIPT_PATH + "generate_gcd_snapshot.py"
cmd = [script, "-m" , str(MJD), "-l", logfilename, "-s", options.SEASON]
print(cmd)
err_code = subprocess.call(cmd)
print("...done pulling GCD file. %s" % "SUCCESS" if err_code == 0 else "FAILURE")
if err_code != 0:
    sys.exit(err_code)

################################################################################
# Now correct the GCD file that was just generated                             #
################################################################################
print("Correcting the GCD file...")
script = SCRIPT_PATH + "correct_GCD.py"
err_code = subprocess.call([script,
                    "-i", "./gcd_snapshot.i3.gz", 
                    "-o", out_filename, 
                    "-l", logfilename])
print("...done correcting GCD file. %s" \
      % "SUCCESS" if err_code == 0 else "FAILURE")
if err_code != 0:
    sys.exit(err_code)

################################################################################
# Now test it.                                                                 #
################################################################################
print("Testing the GCD file...")
script = SCRIPT_PATH + "test_GCD.py"
err_code = subprocess.call([script, "-i" , out_filename])                    
print("...done. %s" % "SUCCESS" if err_code == 0 else "FAILURE")
if err_code != 0:
    sys.exit(err_code)

################################################################################
# And stress test it.                                                          #
################################################################################
print("Generating a stress-test sample...")
script = SCRIPT_PATH + "generate_IC86_stress_test_samples.py"
err_code = subprocess.call([script,
                    "-g", out_filename,
                    "-f", "stress_test_samples.i3.gz"])
print("...done. %s" % "SUCCESS" if err_code == 0 else "FAILURE")
if err_code != 0 : sys.exit(err_code)

################################################################################
# And finally validate it.                                                     #
################################################################################
print("Validating the stress-test sample...")
script = SCRIPT_PATH + "validate_stress_test_samples.py"
err_code = subprocess.call([script,
                    "-i", "stress_test_samples.i3.gz"])
print("...done. %s" % "SUCCESS" if err_code == 0 else "FAILURE")
if err_code != 0 : 
    sys.exit(err_code)

################################################################################
# Write the svn info used to generate this file.                               #      
################################################################################
cmd = ["svn","info",expandvars("$I3_SRC")]
svn_info = subprocess.check_output(cmd)
logfile.write(svn_info)

################################################################################
# Write the timestamp to the logfile.                                          #
################################################################################
import datetime
timestamp = str(datetime.datetime.now())
logfile.write("TIMESTAMP : %s\n" % timestamp)

################################################################################
# Write the checksum to the logfile.                                           #
################################################################################
cmd = ["md5sum", out_filename]
check_sum = subprocess.check_output(cmd)
logfile.write("MD5SUM : %s\n" % check_sum)

