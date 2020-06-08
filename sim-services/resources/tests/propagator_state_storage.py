#!/usr/bin/env python

"""
Demonstrate that the results of muon propagation and cascade simulation
can be exactly recovered by storing the random number generator state
and re-playing the simulation for a subset of events with the same
random number sequence.
"""

from icecube import icetray, dataclasses, dataio
from I3Tray import I3Tray
from icecube import MuonGun, phys_services, sim_services, simclasses
from icecube.MuonGun.segments import GenerateBundles
from os.path import expandvars
from os import unlink

#if not compiled with sprng support don't do anything
if not hasattr(phys_services,"I3SPRNGRandomService"):
        exit(0)

gcd = expandvars("$I3_TESTDATA/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz")

def make_propagators():
	"""
	Set up a stable of propagators for all the kinds of particles we're going to see.
	"""
	from icecube.PROPOSAL import I3PropagatorServicePROPOSAL
	from icecube.cmc import I3CascadeMCService
	propagators = sim_services.I3ParticleTypePropagatorServiceMap()
	muprop = I3PropagatorServicePROPOSAL()
	cprop = I3CascadeMCService(phys_services.I3GSLRandomService(1)) # dummy RNG
	# should one also consider taus?
	for pt in 'MuMinus', 'MuPlus':
		propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = muprop
	for pt in 'DeltaE', 'Brems', 'PairProd', 'NuclInt', 'Hadrons', 'EMinus', 'EPlus':
		propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = cprop
	return propagators

def generate(nevents=1, fname='foo.i3'):
	tray = I3Tray()
	
	generator = MuonGun.Floodlight()
	
	# set up a random number generator
	randomService = phys_services.I3SPRNGRandomService(
	    seed = 1,
	    nstreams = 10000,
	    streamnum = 1)
	tray.context['I3RandomService'] = randomService

	def make_particle():
		p = dataclasses.I3Particle()
		p.pos = dataclasses.I3Position(0,0,2e3)
		p.dir = dataclasses.I3Direction(0,0)
		p.time = 0
		p.energy = 10**randomService.Uniform(3, 6)
		p.type = p.DeltaE
		p.location_type = p.InIce
		return p
	make_particle.i = 0

	def make_mctree(frame):
		mctree = dataclasses.I3MCTree()
		primary = make_particle()
		primary.location_type = primary.Anywhere
		primary.type = primary.unknown
		muon = make_particle()
		mctree.add_root(primary)
		mctree.append_child(primary, muon)
		frame['I3MCTree'] = mctree
	
	tray.AddSegment(GenerateBundles, 'BundleGen', Generator=generator,
	    NEvents=nevents, GCDFile=gcd)	
	
	def stash(frame):
		frame['RemadeMCTree'] = dataclasses.I3MCTree(frame['I3MCTree'])
	tray.AddModule(stash, 'copy', Streams=[icetray.I3Frame.DAQ])
	
	tray.AddModule('I3PropagatorModule', 'propagator', PropagatorServices=make_propagators(),
	    RandomService=randomService, RNGStateName="RNGState")
	
	tray.AddModule('I3Writer', 'writer',
	    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
	    filename=fname)
	
	
	tray.Execute(nevents+3)
	
	
def check(fname='foo.i3', fraction=0.1):
	
	from unittest import TestCase
	
	tray = I3Tray()
	
	tray.AddModule('I3Reader', 'reader', filenamelist=[gcd, fname])
	
	# set up a random number generator
	randomService = phys_services.I3SPRNGRandomService(
	    seed = 2,
	    nstreams = 10000,
	    streamnum = 1)
	
	# Throw out a random subset of events. The remainder must be reproducible
	# given the same random number sequence
	def drop_random_events(frame, fraction=fraction):
		return randomService.uniform(0, 1) < (1-fraction)
	tray.AddModule(drop_random_events, 'triggerhappy', Streams=[icetray.I3Frame.DAQ])
	
	tray.Add("Rename", keys=["MMCTrackList", "OldMMCTrackList"])
	
	tray.AddModule('I3PropagatorModule', 'propagator', PropagatorServices=make_propagators(),
	    RandomService=randomService, RNGStateName="RNGState", InputMCTreeName="RemadeMCTree", OutputMCTreeName="RemadeMCTree")
	
	class MCTreeCompare(TestCase):
		"""
		Ensure that every entry in the re-simulated MCTree is completely
		identical to the original one.
		"""
		def setUp(self):
			self.orig_tree = self.frame['I3MCTree']
			self.orig_mmctracks = self.frame['OldMMCTrackList']
			self.new_tree = self.frame['RemadeMCTree']
			self.new_mmctracks = self.frame['MMCTrackList']
			
		def testTotalSize(self):
			self.assertEquals(len(self.orig_tree), len(self.new_tree))
			self.assertEquals(len(self.orig_mmctracks), len(self.new_mmctracks))
		def testParticleContent(self):
			
			for p1, p2 in zip(self.orig_tree, self.new_tree):
				if p1.location_type != p1.InIce:
					continue
				self.assertEquals(p1.energy, p2.energy)
				self.assertEquals(p1.type, p2.type)
				self.assertEquals(p1.time, p2.time)
				self.assertEquals(p1.length, p2.length)
				for d in 'x', 'y', 'z':
					self.assertEquals(getattr(p1.pos, d), getattr(p2.pos, d))
				for d in 'zenith', 'azimuth':
					self.assertEquals(getattr(p1.dir, d), getattr(p2.dir, d))
			for t1, t2 in zip(self.orig_mmctracks, self.new_mmctracks):
				for prop in 'E', 'X', 'Y', 'Z', 'T':
					for location in 'i', 'c', 'f':
						method = 'Get'+prop+location
						self.assertEquals(getattr(t1, method)(), getattr(t2, method)()) 
	
	tray.AddModule(icetray.I3TestModuleFactory(MCTreeCompare), 'testy', Streams=[icetray.I3Frame.DAQ])
	
	
	tray.Execute()
	
	

fname = 'propagator_state_storage.i3'
generate(nevents=int(1e2), fname=fname)
check(fraction=0.9, fname=fname)

unlink(fname)

