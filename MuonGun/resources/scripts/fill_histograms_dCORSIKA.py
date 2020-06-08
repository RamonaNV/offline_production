#!/usr/bin/env python

from __future__ import print_function

from icecube import icetray, dataclasses, dataio
from I3Tray import I3Tray

import sys
infiles = sys.argv[1:-1]
outfile = sys.argv[-1]

tray = I3Tray()

tray.AddModule('I3Reader', 'reader', filenamelist=infiles)

import numpy
from icecube.icetray import I3Units

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

import dashi, numpy, tables
from icecube.dataclasses import I3Constants
from icecube.phys_services import I3Calculator
from icecube import sim_services, simclasses, MuonGun
cylinder_intersection = MuonGun.I3MUSICService.cylinder_intersection

class CylinderWeighter(object):
	"""
	In VOLUMECORR mode, CORSIKA generates showers with a zenith distribution proportional
	to the projected area of the cylindrical sampling surface. To convert counts of muons
	detected at the sampling surface back to fluxes, we have to weight each entry
	by 1/(dA/dcos(theta)), the differential projected area of the sampling surface
	into a plane perpendicular to the shower axis.
	"""
	def __init__(self, radius, height, timescale, depthbins):
		from numpy import pi, sqrt
		chi = 2*height/(pi*radius)
		# convert depths (in kmwe) to z coordinates
		z = I3Constants.SurfaceElev - I3Constants.OriginElev - (depthbins*I3Units.km/(0.917*1.005))
		def diffarea(ct, r, l):
			"""
			differential projected area of an upright cylinder: dA_perp d\Omega/(dcos(theta))
			"""
			return 2*pi**2*r*(r*ct + (2*l/pi)*sqrt(1-ct**2))
		def bandarea(ct, r, l):
			return diffarea(ct, r, l) - diffarea(ct, r, 0)
		def intarea(r, l):
			"""
			projected area of an upright cylinder, integrated over the upper half-sphere
			"""
			return (pi**2*r*(r+l))
		class weighter(object):
			def __init__(self, f, timescale, r, l):
				self.f = f
				self.timescale = timescale
				self.r = r
				self.l = l
			def __call__(self, ct):
				return self.timescale*self.f(ct, self.r, self.l)
		self.weights = []
		for zhi, zlo in zip(z[:-1], z[1:]):
			print(zlo, zhi, end=' ')
			if zlo > height/2. or zhi < -height/2.:
				# outside the cylinder
				weight = None
			elif zhi > height/2. and zlo <= height/2.:
				# this layer includes the top of the cylinder
				sideheight = (height/2.-zlo)
				print(sideheight, end=' ')
				weight = weighter(diffarea, timescale, radius, sideheight)
			else:
				# a side layer has no top surface
				if zlo <= -height/2.:
					# XXX HACK: ucr's target cylinder is displaced slightly w.r.t 
					# the IceCube coordinate system. Adjust the effective area of
					# the bottom-most slice to compensate.
					sideheight = zhi+height/2. - 5.
				else:
					sideheight = zhi-zlo
				weight = weighter(bandarea, timescale, radius, sideheight)
				
			print('')
			self.weights.append(weight)
		
	def __call__(self, depthidx, zenith):
		wt = 1./self.weights[depthidx](numpy.cos(zenith))
		#print depthidx, self.weights[depthidx][0]
		return wt

class Filla(icetray.I3Module):
	def __init__(self, ctx):
		icetray.I3Module.__init__(self, ctx)
		self.AddOutBox("OutBox")
	
	def Configure(self):
		from collections import defaultdict
		
		depthbins = numpy.linspace(1.0, 5.0, 15)
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
		primary = frame['I3MCTree'].primaries[0]
		
		# Bail if no muon made it to the sampling surface
		mmctracks = frame['MMCTrackList']
		if len(mmctracks) == 0:
			self.PushFrame(frame)
			return
		
		if self.weighter is None:
			# generated by I3CORSIKAWeightModule
			wm = frame['CorsikaWeightMap']
			print(dict(wm))
			if wm['Weight'] != 1.0:
				raise ValueError("Can't deal with weighted DCorsika")
			self.timescale = wm['TimeScale']
			self.cylinder_radius = wm['CylinderRadius']
			self.cylinder_height = wm['CylinderLength']
			self.area = wm['AreaSum']
			self.weighter = CylinderWeighter(self.cylinder_radius, self.cylinder_height, self.timescale, self.depthbins)
		
		zenith = primary.dir.zenith
		zi = max(numpy.searchsorted(self.zenbins, zenith) - 1, 0)
		zif = max(numpy.searchsorted(self.zenbins_fine, zenith) - 1, 0)
		
		multiplicity=self.multiplicity_slices[zif]
		radius=self.radius_slices[zi]
		energy=self.energy_slices[zi]
		
		# Find the point where the shower axis intersects the sampling surface
		intersection = cylinder_intersection(primary.pos, primary.dir, self.cylinder_height, self.cylinder_radius)
		z = primary.pos.z + intersection.first*primary.dir.z
		kmwe = ((I3Constants.SurfaceElev - I3Constants.OriginElev) - z)*(0.917*1.005)/I3Units.km
		#print 'z = %.1f = %.3f kmwe (d = %.1f)' % (z, kmwe, intersection.first)
		di = max(numpy.searchsorted(self.depthbins, kmwe) - 1, 0)
		
		# We know which depth bin we hit, and the projected area of the
		# associated cylinder segment. Calculate a weight that turns counts
		# back into fluxes.
		#print di, numpy.cos(zenith), weight
		
		self.primary.fill_single((zenith, primary.energy), 1./(self.timescale*self.area))
		

		
		weight = self.weighter(di, zenith)
		
		mult = max(len(mmctracks), 1)
		
		multiplicity[di].fill_single(mult, weight)
		
		values = numpy.asarray([(mult, I3Calculator.closest_approach_distance(primary, dataclasses.I3Position(p.GetXi(), p.GetYi(), p.GetZi())), p.GetEi()) for p in mmctracks])
		radius[di].fill(values[:,:2], weight)
		energy[di].fill(values, weight)
		
		self.nevents += 1
		if self.nevents % 1000 == 0:
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

