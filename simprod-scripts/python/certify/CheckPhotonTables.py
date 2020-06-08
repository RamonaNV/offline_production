#!/usr/bin/env python

import os,sys
from os.path import expandvars, exists, join, basename, isdir
import glob
import pickle
import random
import logging
from optparse import OptionParser


usage = "usage: %prog [options]"
parser = OptionParser(usage)

# Photon tables top directory for us, will probably change for you
parser.add_option("-T", "--TableDir", default='/data/sim/sim-new/PhotonTables_SPICEMie_IC79_IC86.2011',
                  dest="TableDir", help="Photon tables top directory")


parser.add_option("-P", "--PickleChecksums", default='/data/sim/sim-new/downloads/PhotonTablesChkSums.dat',
                  dest="ChkSumFile", help="Pickled file containing photon tables and their check sums")

parser.add_option("-d", "--debug", action='store_true', default=False,
                  dest="DEBUG", help="prints out file info when in debug mode, default is false")


(options,args) = parser.parse_args()
if len(args) != 0:
    message = "Got undefined options:"
    for a in args:
        message += a
        message += " "
    parser.error(message)
    
PhotonsTableDir = options.TableDir
if not os.path.isdir(PhotonsTableDir): raise Exception("input Photon Tables directory does not exist")
	
CFile = options.ChkSumFile
if not os.path.isfile(CFile) : raise Exception("input file containing checksums does not exist")

#Load pickled file containing random selection of files and their chksums
Dirs = pickle.load(open(CFile, 'rb'))

debug = options.DEBUG

# List of subdirectories to be checked
PhotorecCscd_ = join(PhotonsTableDir,"AHA07v1ice/L1_showers/")
PhotorecCscdDrivers_ = join(PhotonsTableDir,"AHA07v1ice/driverfiles/")

PhotorecMu_ = join(PhotonsTableDir,"AHA07v2h2/L1_showers/")
PhotorecMuDrivers_ = join(PhotonsTableDir,"AHA07v2h2/driverfiles/")

PhotorecSpice1FiniteReco_ = join(PhotonsTableDir,"SPICE1/L1_showers/")
PhotorecSpice1Mu_ = join(PhotonsTableDir,"SPICE1/L2_muons/")
PhotorecSpice1Drivers_ = join(PhotonsTableDir,"SPICE1/driverfiles/")

PhotonSpiceMieL1_ = join(PhotonsTableDir,"SPICEMie/L1_showers/")
PhotonSpiceMieL2_ = join(PhotonsTableDir,"SPICEMie/L2_muons/")
PhotonSpiceMieDrivers_ = join(PhotonsTableDir,"SPICEMie/driverfiles/")


def md5sum(filename,buffersize=16384):
    """Return md5 digest of file"""
    filed = open(filename)
    try:       import hashlib
    except ImportError:
       import md5
       digest = md5.new()
    else:
       digest = hashlib.md5()

    buffer = filed.read(buffersize)
    while buffer:
        digest.update(buffer)
        buffer = filed.read(buffersize)
    filed.close()
    return digest.hexdigest()

Groups_ = [PhotorecCscd_,PhotorecCscdDrivers_,PhotorecMu_,PhotorecMuDrivers_,PhotorecSpice1FiniteReco_,PhotorecSpice1Mu_,
	   PhotorecSpice1Drivers_,PhotonSpiceMieL1_,PhotonSpiceMieL2_,PhotonSpiceMieDrivers_]
GroupNames = ["PhotorecCscd","PhotorecCscdDrivers","PhotorecMu","PhotorecMuDrivers","PhotorecSpice1FiniteReco",
	      "PhotorecSpice1Mu","PhotorecSpice1Drivers","PhotonSpiceMieL1","PhotonSpiceMieL2","PhotonSpiceMieDrivers"]

#Groups_ = [PhotonSpiceMieL1_]
#GroupNames = ["PhotonSpiceMieL1"]


def DoChecks(Groups_,GroupNames,debug):
    Check = 1
    for GroupName in GroupNames:
	logging.debug("**** Now checking for tables in group '%s'  *****" % GroupName)
	tmpDict = Dirs[GroupName]
	SFiles = random.sample(list(tmpDict.keys()),min(2,len(tmpDict)))
	for f in SFiles:
	    if os.path.isfile(join(Groups_[GroupNames.index(GroupName)],f)) :
		if not tmpDict[f] == md5sum(join(Groups_[GroupNames.index(GroupName)],f)):    
		    Check = 0
		    logging.debug("file %s failed chksum"%join(Groups_[GroupNames.index(GroupName)],f))
	    else:
		Check = 0
		logging.debug("file %s does not exist"%join(Groups_[GroupNames.index(GroupName)],f)) 
    return Check
	
    
Check = DoChecks(Groups_,GroupNames,debug)

if not Check:
    logging.warning("""The check on at least 1 file failed, run in "debug" mode to see which file(s) are affected""")
    
logging.debug(Check)
