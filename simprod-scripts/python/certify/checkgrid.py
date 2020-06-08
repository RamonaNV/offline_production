#!/usr/bin/env python
################################ Tray 0 ##############################
#   on 2011-07-29
#
######################################################################
import commands,os,sys,string
import glob
import time
import random
import logging

from threading import Thread
from os.path import expandvars

from I3Tray import *
from iceprod.modules import ipmodule

class CertifyGrid(ipmodule.IPBaseClass):
    """
    Wrapper class that copy files to storage units with GridFTP
    """
    def __init__(self):
        ipmodule.IPBaseClass.__init__(self)

        self.AddParameter('photontablesdir','main dir containing photon and photorec tables sub dirs','$system(photontablesdir)')
        self.AddParameter('lib_url','downloads directory','$system(lib_url)')
        self.AddParameter('ptversion','photon tables version','$steering(ptversion)')


    def GetCheckSum(self,FileArg):
        import md5
        try:
            FileArgBin = open(FileArg, "rb")
            FileArgBin1 = FileArgBin.read()
            MCipher = md5.new()
            MCipher.update(FileArgBin1)
            CipherHex = MCipher.hexdigest()

            return CipherHex
        except:
            return 0


    def DoChecks(self,DirName,Dict):
        import random
        SFiles = list(Dict.keys())
        for f in SFiles:
            if os.path.isfile(DirName+"/"+f) :
                if not Dict[f] == self.GetCheckSum(DirName+"/"+f):
                    return 0
            else:
                return 0
        return 1


    def Execute(self,stats):
        import datetime
        STime = datetime.datetime.now()
        
        start_ = time.clock()
        
        if not ipmodule.IPBaseClass.Execute(self,stats): return 0

        import icecube.icetray
        import cPickle
        import random
        import md5

        PhotonTablesChkSums = cPickle.load(open('PhotonTablesChkSums.dat', 'rb'))

        TopDir = expandvars("$TMPTOP")
        syschkout = os.path.join(TopDir,"syschk.out")
        outfile = open(syschkout,"w")

        LibUrl = self.GetParameter('lib_url')
        LibUrl = self.parser.parse(LibUrl)
        logging.info("LibUrl is %s "%LibUrl)
        outfile.write("LibUrl=%s\n"%LibUrl)

        ptversion = self.GetParameter('ptversion')
        ptversion = self.parser.parse(ptversion)
        logging.info("ptversion %s "%ptversion)
        outfile.write("ptversion=%s\n"%ptversion)

        PhotonTablesDir = self.GetParameter('photontablesdir')
        PhotonTablesDir = self.parser.parse(PhotonTablesDir)
        logging.info("PhotonTablesDir is %s "%PhotonTablesDir)
        outfile.write("PhotonTablesDir=%s\n"%PhotonTablesDir)

        if os.path.exists(PhotonTablesDir+"/tables"): stats["HasPhotonTablesDir"] = 1

        stats["HasPhotorecCscd_AHA07v1ice"] = 0
        stats["HasPhotorecCscdDriverDir"] = 0
        stats["PhotorecCscdChkSum"] = 0
        stats["PhotorecCscdDriverChkSum"] = 0
        if os.path.exists(PhotonTablesDir+"/AHA07v1ice/L1_showers") :
            stats["HasPhotorecCscd_AHA07v1ice"] = 1
            stats["PhotorecCscdChkSum"] = self.DoChecks(PhotonTablesDir+"/AHA07v1ice/L1_showers",PhotonTablesChkSums["PhotorecCscd"])
        if os.path.exists(PhotonTablesDir+"/AHA07v1ice/driverfiles") :
            stats["HasPhotorecCscdDriverDir"] = 1
            stats["PhotorecCscdDriverChkSum"] = self.DoChecks(PhotonTablesDir+"/AHA07v1ice/driverfiles",PhotonTablesChkSums["PhotorecCscdDrivers"])

        stats["HasPhotorecMu_AHA07v2h2"] = 0
        stats["HasPhotorecMuDriverDir"] = 0
        stats["PhotorecMuChkSum"] = 0
        stats["PhotorecMuDriverChkSum"] = 0
        if os.path.exists(PhotonTablesDir+"/AHA07v2h2/L1_showers") :
            stats["HasPhotorecMu_AHA07v2h2"] = 1
            stats["PhotorecMuChkSum"] = self.DoChecks(PhotonTablesDir+"/AHA07v2h2/L1_showers",PhotonTablesChkSums["PhotorecMu"])
        if os.path.exists(PhotonTablesDir+"/AHA07v2h2/driverfiles") :
            stats["HasPhotorecMuDriverDir"] = 1
            stats["PhotorecMuDriverChkSum"] = self.DoChecks(PhotonTablesDir+"/AHA07v2h2/driverfiles",PhotonTablesChkSums["PhotorecMuDrivers"])

        stats["HasPhotorecFinitereco_SPICE1_L1"] = 0
        stats["HasPhotorecMu_SPICE1_L2"] = 0
        stats["HasPhotorecSPICE1DriverDir"] = 0
        stats["PhotorecFinitereco_SPICE1_L1ChkSum"] = 0
        stats["PhotorecMu_SPICE1_L2ChkSum"] = 0
        stats["PhotorecSPICE1DriverChkSum"] = 0
        if os.path.exists(PhotonTablesDir+"/SPICE1/L1_showers/") :
            stats["HasPhotorecFinitereco_SPICE1_L1"] = 1
            stats["PhotorecFinitereco_SPICE1_L1ChkSum"] = self.DoChecks(PhotonTablesDir+"/SPICE1/L1_showers",PhotonTablesChkSums["PhotorecSpice1FiniteReco"])
        if os.path.exists(PhotonTablesDir+"/SPICE1/L2_muons/") :
            stats["HasPhotorecMu_SPICE1_L2"] = 1
            stats["PhotorecMu_SPICE1_L2ChkSum"] = self.DoChecks(PhotonTablesDir+"/SPICE1/L2_muons",PhotonTablesChkSums["PhotorecSpice1Mu"])
        if os.path.exists(PhotonTablesDir+"/SPICE1/driverfiles/") :
            stats["HasPhotorecSPICE1DriverDir"] = 1
            stats["PhotorecSPICE1DriverChkSum"] = self.DoChecks(PhotonTablesDir+"/SPICE1/driverfiles",PhotonTablesChkSums["PhotorecSpice1Drivers"])

        stats["HasPhotonTables_L1_showers"] = 0
        stats["HasPhotonTables_L2_muons"] = 0
        stats["HasPhotonDriverDir"] = 0
        stats["PhotonTables_L1_showersChkSum"] = 0
        stats["PhotonTables_L2_muonsChkSum"] = 0
        stats["PhotonDriverChkSum"] = 0
        if os.path.exists(PhotonTablesDir+"/"+ptversion+"/L1_showers") :
            stats["HasPhotonTables_L1_showers"] = 1
            stats["PhotonTables_L1_showersChkSum"] = self.DoChecks(PhotonTablesDir+"/"+ptversion+"/L1_showers/",PhotonTablesChkSums["PhotonSpiceMieL1"])
        if os.path.exists(PhotonTablesDir+"/"+ptversion+"/L2_muons") :
            stats["HasPhotonTables_L2_muons"] = 1
            stats["PhotonTables_L2_muonsChkSum"] = self.DoChecks(PhotonTablesDir+"/"+ptversion+"/L2_muons/",PhotonTablesChkSums["PhotonSpiceMieL2"])
        if os.path.exists(PhotonTablesDir+"/"+ptversion+"/driverfiles") :
            stats["HasPhotonDriverDir"] = 1
            stats["PhotonDriverChkSum"] = self.DoChecks(PhotonTablesDir+"/"+ptversion+"/driverfiles/",PhotonTablesChkSums["PhotonSpiceMieDrivers"])

        stats["PhotonTablesOk"] = stats["PhotorecCscdChkSum"] * stats["PhotorecCscdDriverChkSum"] * stats["PhotorecMuChkSum"] * stats["PhotorecMuDriverChkSum"] \
                                * stats["PhotorecFinitereco_SPICE1_L1ChkSum"] * stats["PhotorecMu_SPICE1_L2ChkSum"] * stats["PhotorecSPICE1DriverChkSum"] \
                                * stats["PhotonTables_L1_showersChkSum"] * stats["PhotonTables_L2_muonsChkSum"] * stats["PhotonDriverChkSum"]


        stats["HasPhotonTablesSPICE1_L1_showers"] = 0
        stats["HasPhotonTablesSPICE1_L2_muons"] = 0
        stats["HasPhotonDriverDirSPICE1"] = 0
        stats["PhotonTablesSPICE1_L1_showersChkSum"] = 0
        stats["PhotonTablesSPICE1_L2_muonsChkSum"] = 0
        stats["PhotonDriverChkSumSPICE1"] = 0
        if os.path.exists(PhotonTablesDir+"/SPICE1/L1_showers") :
            stats["HasPhotonTablesSPICE1_L1_showers"] = 1
            stats["PhotonTablesSPICE1_L1_showersChkSum"] = self.DoChecks(PhotonTablesDir+"/SPICE1/L1_showers/",PhotonTablesChkSums["PhotonSpice1L1"])
        if os.path.exists(PhotonTablesDir+"/SPICE1/L2_muons") :
            stats["HasPhotonTablesSPICE1_L2_muons"] = 1
            stats["PhotonTablesSPICE1_L2_muonsChkSum"] = self.DoChecks(PhotonTablesDir+"/SPICE1/L2_muons/",PhotonTablesChkSums["PhotonSpice1L2"])
        if os.path.exists(PhotonTablesDir+"/SPICE1/driverfiles") :
            stats["HasPhotonDriverDirSPICE1"] = 1
            stats["PhotonDriverChkSumSPICE1"] = self.DoChecks(PhotonTablesDir+"/SPICE1/driverfiles/",PhotonTablesChkSums["PhotonSpice1Drivers"])


        stats["PhotonTablesSPICE1_Ok"] = stats["HasPhotonTablesSPICE1_L1_showers"] * stats["HasPhotonTablesSPICE1_L2_muons"] * stats["HasPhotonDriverDirSPICE1"] \
                                        * stats["PhotonTablesSPICE1_L1_showersChkSum"] * stats["PhotonTablesSPICE1_L2_muonsChkSum"] * stats["PhotonDriverChkSumSPICE1"]



        stats["Has_cvmfs"] = 0
        if os.path.exists("/cvmfs/atlas.cern.ch") : stats["Has_cvmfs"] = 1

        os.system("uname -a")
        platform = commands.getoutput("uname -a")
        logging.info("platform information: %s"%platform)
        outfile.write("Platform=%s\n"%str(platform))

        stats["TotalMem(GB)"] = 0
        try:
            MemOutput = (commands.getoutput("cat /proc/meminfo | grep MemTotal")).split()
            MemTotal = float(MemOutput[1])/1024
            stats["TotalMem(GB)"] = float("%.2f"%(MemTotal/1000))
        except:
            pass

        stats["MemPerCPU>=2GB"] = 0
        try:
            NoOfCPUs = commands.getoutput("cat /proc/cpuinfo | grep processor | wc -l")
            MemPerCPU = float(MemTotal/float(NoOfCPUs))
            stats["NoOfCPUs"] = float(NoOfCPUs)
            stats["MemPerCPU(GB)"] = float("%.2f"%(MemPerCPU/1000))
            if MemPerCPU >= 1900.00: stats["MemPerCPU>=2GB"] = 1
        except:
            pass

        stats["CPUSpeed>=2GHz"] = 0
        try:
            CPUOutput = (commands.getoutput("cat /proc/cpuinfo | grep -e 'cpu MHz'")).split()
            CPUModel = (commands.getoutput("cat /proc/cpuinfo | grep -e 'model name'")).split()
            outfile.write("CPUSpeed: %s\n"%CPUOutput)
            outfile.write("CPUModel: %s\n"%CPUModel)
            CPUSpeed = float(CPUOutput[3])
            stats["CPUSpeed(GHz)"] = float("%.2f"%(CPUSpeed/1000))
            if CPUSpeed >= 2000.00: stats["CPUSpeed>=2GHz"] = 1
        except:
            pass

        stats["gccVersion>=3.4.6"]=0
        try:
            gccOutput = (commands.getoutput("gcc --version | grep gcc")).split()
            gccVersion = gccOutput[2]
            if gccVersion >= "3.4.6" : stats["gccVersion>=3.4.6"] = 1
            logging.info("gcc Version is %s"%str(gccVersion))
            outfile.write("gccVersion=%s\n"%gccVersion)
        except:
            pass

        stats["PythonVersion>=2.6"]=0
        stats["PythonVersion"]=0
        try:
            PythonV = sys.version_info
            stats["PythonVersion"]= float(PythonV[0]+float(PythonV[1])/10.0)
            if PythonV[0]>=2 and PythonV[1]>=6 : stats["PythonVersion>=2.6"] = 1
            logging.info("Python Version is %s"%str(PythonV))
            outfile.write("PythonVersion=%s\n"%str(PythonV))

        except:
            pass

        stats["HasJAVA_HOME"] = 0
        stats["JAVA_Version>=1.5.0"] = 0
        try:
            JAVA_HOME = commands.getoutput("echo $JAVA_HOME")
            if os.path.exists(str(JAVA_HOME)): stats["HasJAVA_HOME"] = 1
            if (stats["HasJAVA_HOME"]):
                JAVA_VERSION = (commands.getoutput("$JAVA_HOME/bin/java -version")).split()
                JAVA_VERSION = JAVA_VERSION[2][1:6]
                if JAVA_VERSION >= "1.5.0" : stats["JAVA_Version>=1.5.0"] = 1
                logging.info("Java Version is %s"%str(JAVA_VERSION))
                outfile.write("Java_Version=%s\n"%JAVA_VERSION)

        except:
            pass

        #### check system speed profile ####
        n = 1000
        maxVal = 1000
        m = []
        p = []
        for i in range(n):
            m.append([random.randint(0,maxVal) for el in range(n)])
            p.append([random.randint(0,maxVal) for el in range(n)])

        C = [[0 for i in range(n)] for j in range(n)]
        C = [[sum(m[i][k]*p[k][j] for k in range(len(p))) for j in range(len(p[0]))] for i in range(len(m))]
        ####################################        

        t_d = time.clock() - start_
        outfile.write("Processing time %s"%str(t_d))
        
        stats["ProcessingTime"] = t_d
        
        sys.stdout.flush()
        sys.stderr.flush()
        
        outfile.close()

        return 0

class CopyLogs(ipmodule.IPBaseClass):
    """
    Copy log/outfiles to storage units with GridFTP
    """
    def Execute(self,stats):
        if not ipmodule.IPBaseClass.Execute(self,stats): return 0

        import icecube.icetray
        from iceprod.modules import gsiftp
        from iceprod.core.dataclasses import I3PostTray

        time.sleep(60)

        sys.stdout.flush()
        sys.stderr.flush()

        post    = I3PostTray()
        job     = int(self.parser.parse("$args(procnum)"))
        key     = self.parser.parse("$args(key)")
        dataset = int(self.parser.parse("$args(dataset)"))
        TopDir = expandvars("$TMPTOP")
        

        station = self.parser.parse("$system(gridname)")
        tx_protocol = "gsiftp://gridftp.icecube.wisc.edu"

        if str(station) in ("npx4","npx3","npx2-uwa"): tx_protocol = "file:"


        sys.stdout.flush()
        sys.stderr.flush()
        os.system(""" tar -cf tarred_file.tar %s/*log* %s/*out* %s/*err*"""%(TopDir,TopDir,TopDir))

        post.AddModule("gsiftp.URLCopy","upload_tar")(
        ("source","file:%s" % os.path.join(os.getcwd(),"tarred_file.tar")),
        ("destination","%s/data/sim/scratch/IceSim/SysChkLogs/%s/%s/%s_tarred_file.tar"%(tx_protocol,str(station),str(dataset),str(job))),
        )

        post.Execute()

        return 0
