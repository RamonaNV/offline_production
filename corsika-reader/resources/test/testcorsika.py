#!/usr/bin/env python

from I3Tray import *
from icecube import icetray, dataio, dataclasses, phys_services, corsika_reader
import math

i3_testdata = os.path.expandvars("$I3_TESTDATA/corsika-reader") \
    if os.environ.get('I3_TESTDATA') \
    else os.path.expandvars("$I3_TESTDATA/corsika-reader")
                

tray = I3Tray()
tray.context['I3RandomService'] = phys_services.I3GSLRandomService(42)
tray.Add(corsika_reader.ReadCorsika, 'reader', FilenameList=[i3_testdata + '/DAT010000'], 
	NEvents=1,OverSampling=10, TrimShower=False)
tray.AddModule('CountFrames', 'counter', DAQ=10)

import unittest
class SanityCheck(unittest.TestCase):
	def testKeys(self):
		self.assert_('I3MCTree' in self.frame, "MC tree generated")
		self.assert_('CorsikaWeightMap' in self.frame, "Weights generated")
	def testPrimary(self):
		primary = dataclasses.get_most_energetic_primary(self.frame['I3MCTree'])
		self.assert_(primary.type == dataclasses.I3Particle.PPlus, "Proton primary")
		self.assert_(abs(primary.dir.zenith - 0.1163) < 0.001, "Zenith matches")
		self.assert_(abs(primary.dir.azimuth - 1.94033) < 0.001, "Azimuth matches")
		self.assert_(abs(primary.energy - 71.147*I3Units.TeV) < 1*I3Units.GeV, "Energy matches")
		self.assert_(math.sqrt(primary.pos.x**2 + primary.pos.y**2 + primary.pos.z**2) < 1600, "Inside detector volume")

	def testTree(self):
		primary = dataclasses.get_most_energetic_primary(self.frame['I3MCTree'])

		print('%d particles in tree' % len(self.frame['I3MCTree']))
		self.assert_(len(self.frame['I3MCTree']) == 1350, "1350 particles in tree")
		allowed_particles = [dataclasses.I3Particle.MuMinus, dataclasses.I3Particle.MuPlus, dataclasses.I3Particle.NuE, dataclasses.I3Particle.NuEBar, dataclasses.I3Particle.NuMu, dataclasses.I3Particle.NuMuBar, dataclasses.I3Particle.NuTau, dataclasses.I3Particle.NuTauBar, dataclasses.I3Particle.PPlus]
		for particle in self.frame['I3MCTree']:
			self.assert_(particle.type in allowed_particles, "Particle type is one that should be in the MC Tree")
			self.assert_(particle.energy <= primary.energy, "Energy conservation")

	def testWeights(self):
		print(self.frame['CorsikaWeightMap'])
		wdict = self.frame['CorsikaWeightMap']
		self.assert_(abs(wdict['EnergyPrimaryMax'] - 100*I3Units.TeV) < 1 *I3Units.GeV, "Max simulation energy")
		self.assert_(abs(wdict['EnergyPrimaryMin'] - 46.5*I3Units.TeV) < 1 *I3Units.GeV, "Min simulation energy")
		self.assert_(wdict['PrimarySpectralIndex'] == -1, "Spectral index")
		self.assert_(wdict['NEvents'] == 1, "NEvents")


tray.AddModule(icetray.I3TestModuleFactory(SanityCheck), 'testy',Streams=[icetray.I3Frame.DAQ])

tray.Execute()
