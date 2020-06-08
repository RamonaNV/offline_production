"""
IceProd Module for Simple Muon Production
"""
import os,sys
from os.path import expandvars

from I3Tray import I3Units, I3Tray
from icecube import icetray, dataio, dataclasses

from .. import ipmodule
from .. import segments 
from ..util import ReadI3Summary, WriteI3Summary

class SimpleCascade(ipmodule.ParsingModule):
   '''
   Injects a cascade randomly in a cylinder centered on, 
   but well within the IC86 detector.
   '''
   
   def __init__(self):
      ipmodule.ParsingModule.__init__(self)
      
      self.AddParameter('outputfile','Output filename','')
      self.AddParameter("seed","RNG seed",0)
      self.AddParameter('nevents','Number of events',0)
      self.AddParameter('FromEnergy','Minimum energy',1.*I3Units.TeV)
      self.AddParameter('ToEnergy','Maximum energy',10.*I3Units.PeV)
      self.AddParameter('HistogramFilename', 'Histogram filename.', None)
      
   def Execute(self,stats):
      if not ipmodule.ParsingModule.Execute(self,stats):
         return 0

      import random
      from math import pi, sqrt, sin, cos

      random.seed(self.seed)
       
      # Instantiate a tray 
      tray = I3Tray()
            
      tray.Add("I3InfiniteSource",
               Stream = icetray.I3Frame.DAQ)

      def Generator(frame, FromEnergy = 1*I3Units.TeV, ToEnergy = 1*I3Units.TeV):
         p = dataclasses.I3Particle()
         p.energy = random.uniform(FromEnergy, ToEnergy)
         p.pos = dataclasses.I3Position(0,0,0)

         # sample on a cylinder
         theta = random.uniform(0., 2*pi)
         r = sqrt(random.uniform(0, 300*I3Units.m * 300*I3Units.m))

         x = r * sin(theta)
         y = r * cos(theta)
         z = random.uniform(-300*I3Units.m, 300*I3Units.m)        
         
         zenith = random.uniform(0., pi)
         azimuth = random.uniform(0., 2*pi)
         p.dir = dataclasses.I3Direction(zenith, azimuth)
         p.length = 500 * I3Units.m
         p.type = dataclasses.I3Particle.ParticleType.EMinus
         p.location_type = dataclasses.I3Particle.LocationType.InIce
         p.time = 0. * I3Units.ns
         
         tree = dataclasses.I3MCTree()
         tree.add_primary(p)
                       
         frame["I3MCTree"] = tree
         
      tray.Add(Generator,
               FromEnergy = self.fromenergy,
               ToEnergy = self.toenergy,
               Streams = [icetray.I3Frame.DAQ]
      )


      if self.histogramfilename:         
         from icecube.production_histograms import ProductionHistogramModule
         from icecube.production_histograms.histogram_modules.simulation.mctree_primary import I3MCTreePrimaryModule
         from icecube.production_histograms.histogram_modules.simulation.mctree import I3MCTreeModule
        
         tray.AddModule(ProductionHistogramModule, 
                        Histograms = [I3MCTreePrimaryModule, I3MCTreeModule],
                        OutputFilename = self.histogramfilename)
      
      tray.Add("I3Writer", 
               filename = self.outputfile,
               Streams = [icetray.I3Frame.TrayInfo,
                          icetray.I3Frame.DAQ,
                          icetray.I3Frame.Stream('S'),
                          icetray.I3Frame.Stream('M')])

      # Execute the Tray
      print(tray)
      tray.Execute(self.nevents)
      
      # Free memory
      del tray
      return 0


