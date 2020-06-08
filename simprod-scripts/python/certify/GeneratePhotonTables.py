
#!/usr/bin/env python
#

import os
from os.path import expandvars, exists, join, basename, isdir
import sys
import md5
import glob
import pickle
import random
import logging

PhotonsTableDir = "/data/sim/sim-new/PhotonTables_SPICEMie_IC79_IC86.2011"

PhotonsTableSimDir = "/data/sim/sim-new/PhotonTablesProduction"

PhotorecCscd_ = join(PhotonsTableDir,"AHA07v1ice/L1_showers/")
PhotorecCscdDrivers_ = join(PhotonsTableDir,"AHA07v1ice/driverfiles/")
PhotorecCscd = {}
PhotorecCscdDrivers = {}

PhotorecMu_ = join(PhotonsTableDir,"AHA07v2h2/L1_showers/")
PhotorecMuDrivers_ = join(PhotonsTableDir,"AHA07v2h2/driverfiles/")
PhotorecMu = {}
PhotorecMuDrivers = {}

PhotorecSpice1FiniteReco_ = join(PhotonsTableDir,"SPICE1/L1_showers/")
PhotorecSpice1Mu_ = join(PhotonsTableDir,"SPICE1/L2_muons/")
PhotorecSpice1Drivers_ = join(PhotonsTableDir,"SPICE1/driverfiles/")
PhotorecSpice1FiniteReco = {}
PhotorecSpice1Mu = {}
PhotorecSpice1Drivers = {}

PhotonSpice1L1_ = join(PhotonsTableSimDir,"SPICE1/L1_showers/")
PhotonSpice1L2_ = join(PhotonsTableSimDir,"SPICE1/L2_muons/")
PhotonSpice1Drivers_ = join(PhotonsTableSimDir,"SPICE1/driverfiles/")
PhotonSpice1L1 = {}
PhotonSpice1L2 = {}
PhotonSpice1Drivers = {}

PhotonSpiceMieL1_ = join(PhotonsTableDir,"SPICEMie/L1_showers/")
PhotonSpiceMieL2_ = join(PhotonsTableDir,"SPICEMie/L2_muons/")
PhotonSpiceMieDrivers_ = join(PhotonsTableDir,"SPICEMie/driverfiles/")
PhotonSpiceMieL1 = {}
PhotonSpiceMieL2 = {}
PhotonSpiceMieDrivers = {}


def GetCheckSum(FileArg):
	try:
		FileArgBin = open(FileArg, "rb")
		FileArgBin1 = FileArgBin.read()
		MCipher = md5.new()
		MCipher.update(FileArgBin1)
		CipherHex = MCipher.hexdigest()
		
		return CipherHex
	except:
		return 0

Groups_ = [PhotorecCscd_,PhotorecCscdDrivers_,PhotorecMu_,PhotorecMuDrivers_,PhotorecSpice1FiniteReco_,PhotorecSpice1Mu_,PhotorecSpice1Drivers_,PhotonSpice1L1_,PhotonSpice1L2_,PhotonSpice1Drivers_,PhotonSpiceMieL1_,PhotonSpiceMieL2_,PhotonSpiceMieDrivers_]
Groups = [PhotorecCscd,PhotorecCscdDrivers,PhotorecMu,PhotorecMuDrivers,PhotorecSpice1FiniteReco,PhotorecSpice1Mu,PhotorecSpice1Drivers,PhotonSpice1L1,PhotonSpice1L2,PhotonSpice1Drivers,PhotonSpiceMieL1,PhotonSpiceMieL2,PhotonSpiceMieDrivers]
GroupNames = ["PhotorecCscd","PhotorecCscdDrivers","PhotorecMu","PhotorecMuDrivers","PhotorecSpice1FiniteReco","PhotorecSpice1Mu","PhotorecSpice1Drivers","PhotonSpice1L1","PhotonSpice1L2","PhotonSpice1Drivers","PhotonSpiceMieL1","PhotonSpiceMieL2","PhotonSpiceMieDrivers"]


CheckSumDicts = {}

for Group_ in Groups_ :
    
    logging.info("processing files from dir: %s"%Group_)
    
    Files = glob.glob(Group_+"*")
    Files = [f for f in glob.glob(Group_+"*") if "md5sum" not in f ]
    Files.sort()
    
    #Files = Files[0:2]
    
    #SFiles = random.sample(Files,min(50,len(Files)))
    SFiles = Files
    
    for SFile in SFiles:
	(Groups[Groups_.index(Group_)])[os.path.basename(SFile)] = GetCheckSum(SFile)

    CheckSumDicts[GroupNames[Groups_.index(Group_)]] = Groups[Groups_.index(Group_)]

FileOutput = open('PhotonTablesChkSums.dat', 'wb')
pickle.dump(CheckSumDicts, FileOutput)
FileOutput.close()

