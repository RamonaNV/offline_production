#!/usr/bin/env python

from icecube import icetray, dataclasses, dataio
from I3Tray import I3Tray

import sys
infiles = sys.argv[1:-1]
outfile = sys.argv[-1]

tray = I3Tray()

tray.AddModule('I3Reader', 'reader', filenamelist=infiles)

import numpy
from icecube.icetray import I3Units

class HoerandelWeight(object):
	gamma = numpy.array([2.71, 2.64, 2.54, 2.75, 2.95, 2.66, 2.72, 2.68, 2.69, 2.64, 2.66, 2.64, 2.66, 2.75, 2.69, 2.55, 2.68, 2.64, 2.65, 2.7, 2.64, 2.61, 2.63, 2.67, 2.46, 2.59])
	flux = numpy.array([0.0873, 0.0571, 0.00208, 0.000474, 0.000895, 0.0106, 0.00235, 0.0157, 0.000328, 0.0046, 0.000754, 0.00801, 0.00115, 0.00796, 0.00027, 0.00229, 0.000294, 0.000836, 0.000536, 0.00147, 0.000304, 0.00113, 0.000631, 0.00136, 0.00135, 0.0204])
	flux *= (I3Units.TeV/I3Units.GeV)**(gamma-1) # unit conversion
	def __init__(self, generation_prob, z=1, nevents=1e6*1):
		self.generation_prob = generation_prob
		zi = int(z)-1
		self.z = z
		self.norm = self.flux[zi]
		self.gamma = self.gamma[zi]
	
	@staticmethod	
	def fluxdiff(e, z, gamma, delta_gamma=2.1, eps_cutoff=1.9, E_knee=4.49*I3Units.PeV):
		"""
		Differential (unnormalized) Hoerandel flux
		"""
		return e**(-gamma)*(1+(e/(E_knee*z))**eps_cutoff)**(-delta_gamma/eps_cutoff)
		
	def __call__(self, E):
		return self.norm*self.fluxdiff(E, self.z, self.gamma)/self.generation_prob(E)

class HoerandelWeight5(HoerandelWeight):
	"""
	Hoerandel with only 5 components, after Becherini et al. (also the same as Arne Schoenwald's version)
	"""
	gamma = numpy.nan*numpy.zeros(26)
	flux  = numpy.nan*numpy.zeros(26)
	gamma[0]  = 2.71 # H
	gamma[1]  = 2.64 # He
	gamma[6]  = 2.58 # N
	gamma[12] = 2.67 # Al
	gamma[25] = 2.58 # Fe
	flux[0]  = 8.73e-2 # H
	flux[1]  = 5.71e-2 # He
	flux[6]  = 3.24e-2 # N
	flux[12] = 3.16e-2 # Al
	flux[25] = 2.18e-2 # Fe
	flux *= (I3Units.TeV/I3Units.GeV)**(gamma-1) # unit conversion
	def __init__(self, *args, **kwargs):
		super(HoerandelWeight5, self).__init__(*args, **kwargs)
		if numpy.isnan(self.norm):
			raise ValueError("I can't account for nuclei with charge %d" % self.z)

import dashi
class buffering_histogram(dashi.histogram.histogram):
	"""
	A histogram that consolidates multiple calls to fill() into more efficient blocks
	"""
	maxbuf = 65535
	def __init__(self, *args, **kwargs):
		super(buffering_histogram, self).__init__(*args, **kwargs)
		self._bh_sample = numpy.empty((self.maxbuf, self.ndim))
		self._bh_weights = numpy.empty(self.maxbuf)
		self._bh_pos = 0
	def flush(self):
		super(buffering_histogram, self).fill(self._bh_sample[:self._bh_pos], self._bh_weights[:self._bh_pos])
		self._bh_pos = 0
	def fill_single(self, values, weight):
		"""
		Add a single entry with the given weight
		"""
		self._bh_sample[self._bh_pos] = values
		self._bh_weights[self._bh_pos] = weight
		self._bh_pos += 1
		if self._bh_pos == self.maxbuf:
			self.flush()
	def fill(self, sample, weight):
		"""
		Add multiple entries, each with 1/N of the given weight
		"""
		n = sample.shape[0]
		# For large entries, kick straight to vanilla fill()
		if n > self.maxbuf:
			super(buffering_histogram, self).fill(sample, numpy.ones(n)*(weight/n))
			return 
		elif self._bh_pos + n > self.maxbuf:
			self.flush()
		self._bh_sample[self._bh_pos:self._bh_pos+n] = sample
		self._bh_weights[self._bh_pos:self._bh_pos+n] = weight/n
		self._bh_pos += n
		if self._bh_pos == self.maxbuf:
			self.flush()

# here we have heterogenous generation spectra, and have to be very
# careful about the number of files generated for with each setting
from icecube.weighting.weighting import PowerLaw
spectra = {
	dataclasses.I3Particle.PPlus :  (1000 + 64)*PowerLaw(1e6, 5e2, 1e10, -2) + 1000*PowerLaw(1e6, 1e4, 1e10, -2) + 521*PowerLaw(1e3, 1e4, 1e10, -1),
	dataclasses.I3Particle.He4Nucleus : 64*PowerLaw(1e6, 2e3, 1e10, -2) + 731*PowerLaw(1e6, 4e4, 1e10, -2) + 525*PowerLaw(1e3, 1e4, 1e10, -1),
	dataclasses.I3Particle.N14Nucleus : 64*PowerLaw(1e6, 7e3, 1e10, -2) + 350*PowerLaw(1e3, 1e4, 1e10, -1),
	dataclasses.I3Particle.Al27Nucleus : 64*PowerLaw(1e6, 13.5e3, 1e10, -2) + 524*PowerLaw(1e3, 1e4, 1e10, -1),
	dataclasses.I3Particle.Fe56Nucleus : 64*PowerLaw(1e6, 3e4, 1e10, -2) + 520*PowerLaw(1e3, 1e4, 1e10, -1),
}

spectra = {
	dataclasses.I3Particle.PPlus :  64*PowerLaw(1e6, 5e2, 1e10, -2),
}

import dashi, numpy, tables
from icecube import sim_services, MuonGun
class Filla(icetray.I3Module):
	def __init__(self, ctx):
		icetray.I3Module.__init__(self, ctx)
		self.AddOutBox("OutBox")
	
	def Configure(self):
		from collections import defaultdict
		
		depthbins = numpy.linspace(1.0, 5.0, 9)
		depthbins -= numpy.diff(depthbins)[0]/2.
		zenbins = numpy.arccos(numpy.linspace(1, 0, 11))
		zenbins_fine = numpy.arccos(numpy.linspace(1, 0, 101))
		multbins = numpy.array([1, 2, 3, 4, 10, 20, 40, 100], dtype=float)
		rbins = numpy.array([0, 5, 10, 15, 25, 45], dtype=float)
		
		self.primary = buffering_histogram(2, (zenbins, numpy.logspace(2, 11, 101)))
		self.multiplicity = dashi.histogram.histogram(3, (zenbins_fine, depthbins, numpy.arange(1, 100)))
		self.radius = dashi.histogram.histogram(4, (zenbins, depthbins, multbins, numpy.linspace(0, numpy.sqrt(250), 101)**2))
		self.energy = dashi.histogram.histogram(5, (zenbins, depthbins, multbins, rbins, numpy.logspace(0, 6, 101)))
		
		self.multiplicity_slices = tuple([tuple([buffering_histogram(1, (numpy.arange(1, 100),)) for j in range(len(depthbins))]) for i in range(len(zenbins_fine))])
		self.radius_slices = tuple([tuple([buffering_histogram(2, (multbins, numpy.linspace(0, numpy.sqrt(250), 101)**2)) for j in range(len(depthbins))]) for i in range(len(zenbins))])
		self.energy_slices = tuple([tuple([buffering_histogram(3, (multbins, rbins, numpy.logspace(0, 6, 101))) for j in range(len(depthbins))]) for i in range(len(zenbins))])
		
		self.depthbins = depthbins
		self.zenbins = zenbins
		self.zenbins_fine = zenbins_fine
		
		self.weighter = None
		
		import os
		if os.path.exists(outfile):
			os.unlink(outfile)
			
		self.nevents = 0
		
	def DAQ(self, frame):
		primary = frame['MCPrimary']
		#if primary.type != primary.PPlus:
		#	return
		# print primary
		if self.weighter is None:
			if 'CorsikaWeightDict' in frame:
				# generated by I3CORSIKAReader
				wm = frame['CorsikaWeightDict']
				if primary.type == primary.PPlus:
					z = 1
				else:
					z = int(primary.type)%100
				genprob = spectra[primary.type]
				self.weighter = HoerandelWeight5(genprob, z)
			elif 'CorsikaWeightMap' in frame:
				# generated by I3CORSIKAWeightModule
				wm = frame['CorsikaWeightMap']
				if wm['Weight'] != 1.0:
					raise ValueError("Can't deal with weighted DCorsika")
				timescale = wm['TimeScale']
				area = wm['AreaSum']
				self.weighter = lambda energy: 1/(timescale*area)
		weight = self.weighter(primary.energy)
		
		zenith = primary.dir.zenith
		zi = max(numpy.searchsorted(self.zenbins, zenith) - 1, 0)
		zif = max(numpy.searchsorted(self.zenbins_fine, zenith) - 1, 0)
		
		self.primary.fill_single((zenith, primary.energy), weight)
		
		multiplicity=self.multiplicity_slices[zif]
		radius=self.radius_slices[zi]
		energy=self.energy_slices[zi]
		
		for di, (depth, tracks) in enumerate(frame['Tracks'].items()):
			kmwe = depth/I3Units.km
			mult = len(tracks)
			values = numpy.asarray([(mult, p.radius, p.energy) for p in tracks])
			
			multiplicity[di].fill_single(mult, weight)
			radius[di].fill(values[:,:2], weight)
			energy[di].fill(values, weight)
		
		self.nevents += 1
		if self.nevents % 10000 == 0:
			print('%d events' % self.nevents)
		
		self.PushFrame(frame)
		
	def Finish(self):
		for i in range(len(self.zenbins_fine)):
			for j in range(len(self.depthbins)):
				self.multiplicity_slices[i][j].flush()
				self.multiplicity._h_bincontent[i+1,j+1,:] += self.multiplicity_slices[i][j]._h_bincontent
				self.multiplicity._h_squaredweights[i+1,j+1,:] += self.multiplicity_slices[i][j]._h_squaredweights
		for i in range(len(self.zenbins)):
			for j in range(len(self.depthbins)):
				self.radius_slices[i][j].flush()
				self.energy_slices[i][j].flush()
				self.radius._h_bincontent[i+1,j+1,:,:] += self.radius_slices[i][j]._h_bincontent
				self.radius._h_squaredweights[i+1,j+1,:,:] += self.radius_slices[i][j]._h_squaredweights
				self.energy._h_bincontent[i+1,j+1,:,:,:] += self.energy_slices[i][j]._h_bincontent
				self.energy._h_squaredweights[i+1,j+1,:,:,:] += self.energy_slices[i][j]._h_squaredweights
		self.primary.flush()
		with tables.openFile(outfile, 'w') as hdf:
			dashi.histsave(self.primary, hdf, '/', 'primary')
			dashi.histsave(self.multiplicity, hdf, '/', 'multiplicity')
			dashi.histsave(self.radius, hdf, '/', 'radius')
			dashi.histsave(self.energy, hdf, '/', 'energy')
		
tray.AddModule(Filla, 'filla')


tray.Execute()

