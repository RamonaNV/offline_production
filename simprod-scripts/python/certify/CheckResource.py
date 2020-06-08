#!/usr/bin/env python
################################ Tray 0 ##############################
#   on 2013-12-19
#
######################################################################
import commands,os,sys,string
from os.path import expandvars
import glob
import random
import time
import cPickle
import socket
import logging
from threading import Thread

import icecube.icetray
from I3Tray import *
from iceprod.modules import ipmodule
from iceprod.modules import gsiftp


class CheckNode(ipmodule.IPBaseClass):
    """
    Wrapper class that copy files to storage units with GridFTP
    """
    def __init__(self):
        ipmodule.IPBaseClass.__init__(self)

        self.AddParameter('photontablesdir','main dir containing photon and photorec tables sub dirs','$system(photontablesdir)')
        self.AddParameter('lib_url','downloads directory','$system(lib_url)')
        self.AddParameter('ptversion','photon tables version','$steering(ptversion)')

        
    def md5sum(self,FileName,buffersize=16384):
        """Return md5 digest of file"""
        filed = open(FileName)
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


    def PhotonicsFileChecks(self,DirName,Dict):
        #import random

        SFiles = list(Dict.keys())
        #SFiles = random.sample(list(Dict.keys()),min(5,len(Dict)))
        for f in SFiles:
            if os.path.isfile(DirName+"/"+f) :
                if not Dict[f] == self.md5sum(DirName+"/"+f):    
                    return 0
            else:
                return 0
        
        return 1
    
    def ProfileCPU(self):

        CPU_Times = []
        
        #for i in range(100):
        for i in range(3):     # reduce to only 3 trials for longer test
            start_ = time.clock()
            #n = 165 # no. that yields about 1 sec on the npx4 benchmarking platform
            n=1000 # much longer test, takes about 261 seconds on fastest npx4 CPU
            maxVal = 1000
            m = []
            p = []
        
            for i in range(n):
                m.append([random.randint(0,maxVal) for el in range(n)])
                p.append([random.randint(0,maxVal) for el in range(n)])

            C = [[0 for i in range(n)] for j in range(n)]
            C = [[sum(m[i][k]*p[k][j] for k in range(len(p))) for j in range(len(p[0]))] for i in range(len(m))]

            #CPU_Times.append(1./(time.clock() - start_))
            CPU_Times.append(261./(time.clock() - start_))

        return sum(CPU_Times)/len(CPU_Times)        #mean
    
    def CopyLogs(self,TopDir,job,dataset,station):
        from iceprod.modules import gsiftp
        from iceprod.core.dataclasses import I3PostTray
        post = I3PostTray()
        
        tx_protocol = "gsiftp://gridftp.icecube.wisc.edu"
        if str(station) in ("npx4","npx3","npx2-uwa","NPX"):
            tx_protocol = "file:"
        
        FileList = glob.glob(os.getcwd()+'/*out')
        FileList.extend(glob.glob(os.getcwd()+'/*log'))
        FileList.extend(glob.glob(os.getcwd()+'/*err'))

        for f in FileList:
            os.system("cp %s %s_%s"%(f,f,str(job)))

        Files2cp = glob.glob(os.getcwd()+'/*%s'%str(job))

        #change 'skua' to str(station) for production
        post.AddModule("gsiftp.URLMultiCopy","upload_logs")(
            ("sourcelist",Files2cp),
            ("destination","%s/data/sim/scratch/IceSim/SysChkLogs/%s/%s/"%(tx_protocol,str(station),str(dataset))),
        )
        
        post.Execute()


    def Execute(self,stats):

        if not ipmodule.IPBaseClass.Execute(self,stats): return 0

        PhotonTablesChkSums = cPickle.load(open('PhotonTablesChkSums.dat', 'rb'))

        TopDir = os.getcwd()
        syschkout = os.path.join(TopDir,"syschk.out")
        outfile = open(syschkout,"w")

        try:
            DomainName = ".".join(socket.gethostbyaddr(socket.gethostname())[0].split(".")[1:])
            outfile.write("DomainName=%s\n"%DomainName)
        except:
            pass
        
        os.system("uname -a")
        platform = commands.getoutput("uname -a")
        logging.debug("platform information: %s"%platform)
        outfile.write("Platform=%s\n"%str(platform))

        LibUrl = self.GetParameter('lib_url')
        LibUrl = self.parser.parse(LibUrl)
        logging.debug("LibUrl is %s "%LibUrl)
        outfile.write("LibUrl=%s\n"%LibUrl)

        ptversion = self.GetParameter('ptversion')
        ptversion = self.parser.parse(ptversion)
        logging.debug("ptversion %s "%ptversion)
        outfile.write("ptversion=%s\n"%ptversion)

        stats["Has_cvmfs"] = 0
        if os.path.isdir("/cvmfs/icecube.opensciencegrid.org"):
            if os.path.isfile("/cvmfs/icecube.opensciencegrid.org/setup.sh"):
                rv_ = os.system("cat /cvmfs/icecube.opensciencegrid.org/setup.sh")
                if not rv_: stats["Has_cvmfs"] = 1   

        PhotonTablesDir = self.GetParameter('photontablesdir')
        PhotonTablesDir = self.parser.parse(PhotonTablesDir)
        logging.debug("PhotonTablesDir is %s "%PhotonTablesDir)
        outfile.write("PhotonTablesDir=%s\n"%PhotonTablesDir)
        
        stats["PhotorecCscdChkSum"] = 0
        stats["PhotorecCscdDriverChkSum"] = 0
        if os.path.exists(PhotonTablesDir+"/AHA07v1ice/L1_showers") :
            stats["PhotorecCscdChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/AHA07v1ice/L1_showers",PhotonTablesChkSums["PhotorecCscd"])
        if os.path.exists(PhotonTablesDir+"/AHA07v1ice/driverfiles") :
            stats["PhotorecCscdDriverChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/AHA07v1ice/driverfiles",PhotonTablesChkSums["PhotorecCscdDrivers"])


        stats["PhotorecMuChkSum"] = 0
        stats["PhotorecMuDriverChkSum"] = 0
        if os.path.exists(PhotonTablesDir+"/AHA07v2h2/L1_showers") :
            stats["PhotorecMuChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/AHA07v2h2/L1_showers",PhotonTablesChkSums["PhotorecMu"])
        if os.path.exists(PhotonTablesDir+"/AHA07v2h2/driverfiles") :
            stats["PhotorecMuDriverChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/AHA07v2h2/driverfiles",PhotonTablesChkSums["PhotorecMuDrivers"])


        stats["PhotorecFinitereco_SPICE1_L1ChkSum"] = 0
        stats["PhotorecMu_SPICE1_L2ChkSum"] = 0
        stats["PhotorecSPICE1DriverChkSum"] = 0
        if os.path.exists(PhotonTablesDir+"/SPICE1/L1_showers/") :
            stats["PhotorecFinitereco_SPICE1_L1ChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/SPICE1/L1_showers",PhotonTablesChkSums["PhotorecSpice1FiniteReco"])
        if os.path.exists(PhotonTablesDir+"/SPICE1/L2_muons/") :
            stats["PhotorecMu_SPICE1_L2ChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/SPICE1/L2_muons",PhotonTablesChkSums["PhotorecSpice1Mu"])
        if os.path.exists(PhotonTablesDir+"/SPICE1/driverfiles/") :
            stats["PhotorecSPICE1DriverChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/SPICE1/driverfiles",PhotonTablesChkSums["PhotorecSpice1Drivers"])


        stats["PhotonTables_L1_showersChkSum"] = 0
        stats["PhotonTables_L2_muonsChkSum"] = 0
        stats["PhotonDriverChkSum"] = 0
        if os.path.exists(PhotonTablesDir+"/"+ptversion+"/L1_showers") :
            stats["PhotonTables_L1_showersChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/"+ptversion+"/L1_showers/",PhotonTablesChkSums["PhotonSpiceMieL1"])
        if os.path.exists(PhotonTablesDir+"/"+ptversion+"/L2_muons") :
            stats["PhotonTables_L2_muonsChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/"+ptversion+"/L2_muons/",PhotonTablesChkSums["PhotonSpiceMieL2"])
        if os.path.exists(PhotonTablesDir+"/"+ptversion+"/driverfiles") :
            stats["PhotonDriverChkSum"] = self.PhotonicsFileChecks(PhotonTablesDir+"/"+ptversion+"/driverfiles/",PhotonTablesChkSums["PhotonSpiceMieDrivers"])

        stats["PhotonTablesOk"] = stats["PhotorecCscdChkSum"] * stats["PhotorecCscdDriverChkSum"] * stats["PhotorecMuChkSum"] * stats["PhotorecMuDriverChkSum"] \
                                * stats["PhotorecFinitereco_SPICE1_L1ChkSum"] * stats["PhotorecMu_SPICE1_L2ChkSum"] * stats["PhotorecSPICE1DriverChkSum"] \
                                * stats["PhotonTables_L1_showersChkSum"] * stats["PhotonTables_L2_muonsChkSum"] * stats["PhotonDriverChkSum"]

    
        try:        
            stats["cpu_speed"] = self.ProfileCPU()
            logging.debug(stats["cpu_speed"])
            outfile.write("CPU Speed %s"%str(stats["cpu_speed"]))
        except:
            pass

        sys.stdout.flush()
        sys.stderr.flush()
        
        outfile.close()
        
        try:
            job = int(self.parser.parse("$args(procnum)"))
            dataset = int(self.parser.parse("$args(dataset)"))
            station = self.parser.parse("$system(gridname)")
            self.CopyLogs(TopDir,job,dataset,station)
        except Exception as err:
            logging.debug("Error: %s\n"%str(err))
        
        logging.debug(stats)
        
        return 0

