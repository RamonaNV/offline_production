import os, re
import logging
from icecube.simprod import ipmodule


class FixCorsikaWeights(ipmodule.ParsingModule):
   """
   Correct weights for misconfigured CORSIKA 5comp sets
   """
   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('nshowers','Number of generated CR showers',1)
        self.AddParameter('procnum','process number',0)
        self.AddParameter('seed','RNG seed',1)
        self.AddParameter('nproc','Number of processes for (RNG)',1)
        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('outputfile','Output filename','corsika.i3')
        self.AddParameter('inputfilelist',"Input filename list",[])
        self.AddParameter('RunCorsika','Run CORSIKA or only generate INPUTS file',False)
        self.AddParameter('mjd','MJD to set the event times',55697)
        self.AddParameter('eslope','CR spectral index (only if ranpri=0)',-2.7)
        self.AddParameter('ranpri','CR spectrum: 0=individual-nuclei, 1=Wiebel-Sooth, 2=Hoerandel, 3=5-component',2)
        self.AddParameter('crtype','CR Particle Type (only if not dcorsika)',0)
        self.AddParameter('pnorm','5-component relative contribution H,He,N,Al,Fe',[10.,5.,3.,2.,1.])
        self.AddParameter('pgam','5-component spectral indices H,He,N,Al,Fe',[2.0,2.0,2.0,2.0,2.0])
        self.AddParameter("RunId","Configure run ID",0)
        self.AddParameter("locut","Enables skew angle cutfoff",1.58)
        self.AddParameter("kcut","minimum neutrino energy required to keep the shower",0)

        self.AddParameter('CORSIKAseed','CORSIKA seed',1)
        self.AddParameter("dslope","Change in spectral index",0)
        self.AddParameter("eprimarymax",'CR max energy',1e5)
        self.AddParameter("eprimarymin","CR min energy",600)
        self.AddParameter("fluxsum","", 0.131475115)
        self.AddParameter("length","",1600)
        self.AddParameter("oversampling","",1)
        self.AddParameter("radius","",800)
        
        self.AddParameter('corsikaVersion','version of corsika to run','v6900')
        self.AddParameter("corsikaName","Corsika binary name","dcorsika")

        self.AddParameter('cthmin','Min theta of injected cosmic rays',0.0)  
        self.AddParameter('cthmax','Max theta of injected cosmic rays',89.99)  
  
        self.AddParameter('ecuts1','hadron min energy (see corsika docs)',273)  
        self.AddParameter('ecuts2','muon min energy (see corsika docs)',273)  
        self.AddParameter('ecuts3','electron min energy (see corsika docs)',0.003)  
        self.AddParameter('ecuts4','photon min energy (see corsika docs)',0.003)  
        self.AddParameter("CutoffType","Sets SPRIC=T (EnergyPerNucleon) or F (EnergyPerParticle) ","EnergyPerNucleon")
        self.AddParameter("UpperCutoffType","Upper cutoff type (defaults to CutoffType)",None)


 
   def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        from I3Tray import I3Tray
        from icecube import icetray,phys_services, dataio, dataclasses
        from icecube.simprod.util import BasicCounter
        from icecube.simprod.generators.weights import Corsika5CompWeightModule, PolygonatoWeightModule
        from icecube.icetray import I3Units
   
        if self.cutofftype == "EnergyPerNucleon" : self.spric = True
        elif self.cutofftype == "EnergyPerParticle" : self.spric = False
        else: raise Exception, "Undefined CutoffType %s" % cutoff_typ
        if not self.uppercutofftype:
            self.uppercutofftype = self.cutofftype
        
        
        # Instantiate a tray 
        tray = I3Tray()

        # Configure IceTray services
        randomService = phys_services.I3SPRNGRandomService(self.seed, self.nproc, self.procnum)
        tray.context["I3RandomService"] = randomService

        tray.AddModule("I3Reader","reader", filenamelist = [self.gcdfile]+self.inputfilelist)

        pnormsum = float(sum(self.pnorm))
        tray.AddModule(Corsika5CompWeightModule,"5compCorsikaweight",
          name                       = "CorsikaWeightMap",
          nevents                    = self.nshowers,
          spric                      = self.spric,
          ThetaMin                   = self.cthmin*I3Units.degree,
          ThetaMax                   = self.cthmax*I3Units.degree,
          cylinderLength             = self.length*I3Units.meter,
          cylinderRadius             = self.radius*I3Units.meter,
          energyprimarymin           = self.eprimarymin*I3Units.GeV,
          energyprimarymax           = self.eprimarymax*I3Units.GeV,

          PrimaryNormalizationH      =  self.pnorm[0]/pnormsum,
          PrimaryNormalizationHe     =  self.pnorm[1]/pnormsum,
          PrimaryNormalizationCNO    =  self.pnorm[2]/pnormsum,
          PrimaryNormalizationMgAlSi =  self.pnorm[3]/pnormsum,
          PrimaryNormalizationFe     =  self.pnorm[4]/pnormsum,

          PrimarySpectralIndexH      = -self.pgam[0],
          PrimarySpectralIndexHe     = -self.pgam[1],
          PrimarySpectralIndexCNO    = -self.pgam[2],
          PrimarySpectralIndexMgAlSi = -self.pgam[3],
          PrimarySpectralIndexFe     = -self.pgam[4],

          OverSamplingFactor         = self.oversampling
        )
        tray.AddModule(PolygonatoWeightModule,"polygonato")

        tray.AddModule(BasicCounter,"count_g", Streams = [icetray.I3Frame.DAQ], 
              name = "Corrected Events", Stats = stats)
       
        tray.AddModule("I3Writer","writer", filename = self.outputfile, streams =[icetray.I3Frame.DAQ] )

        

        # Execute the Tray
        tray.Execute()
        

        del tray
        return 0


if __name__ == "__main__":
   logging.basicConfig()
   rootLogger = logging.getLogger('')
   rootLogger.setLevel(logging.INFO)

   stats = {}
   fixit = FixCorsikaWeights()
   fixit.ExecuteOpts(stats)
   print stats
