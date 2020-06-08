#!/usr/bin/env python

"""
Run MuonGun with a small target surface surrounding DeepCore, and
plot the generated tracks to illustrate the part of the detector
volume that goes un-simulated.
"""

from argparse import ArgumentParser
from os.path import expandvars
parser = ArgumentParser()
parser.add_argument("-g", "--gcd", default=expandvars('$I3_TESTDATA/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'))
parser.add_argument("outfile", help="save plot to file")
args = parser.parse_args()

from icecube import icetray, dataclasses, dataio
from icecube import phys_services, simclasses, MuonGun
from I3Tray import I3Tray
from os.path import expandvars

tray = I3Tray()

tray.context['I3RandomService'] = phys_services.I3GSLRandomService(1337)

from icecube.MuonGun.segments import GenerateBundles
outer = MuonGun.Cylinder(1600, 800)
inner = MuonGun.Cylinder(300, 150, dataclasses.I3Position(0,0,-350))
spectrum = MuonGun.OffsetPowerLaw(5, 1e3, 1e1, 1e4)
model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
generator = MuonGun.EnergyDependentSurfaceInjector(outer, model.flux, spectrum, model.radius,
    MuonGun.ConstantSurfaceScalingFunction(inner))
tray.AddSegment(GenerateBundles, 'BundleGen', NEvents=1000, Generator=generator,
    GCDFile=expandvars('$I3_TESTDATA/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz'))


class Harvest(icetray.I3ConditionalModule):
	def __init__(self, context):
		super(Harvest, self).__init__(context)
		self.AddOutBox("OutBox")
	
	def Configure(self):
		self.tracks = []
	
	def DAQ(self, frame):
		self.tracks.append(dataclasses.get_most_energetic_track(frame['I3MCTree']))
		self.PushFrame(frame)
	
	@staticmethod
	def cylinder_patch(cyl, view='side', **kwargs):
		from matplotlib.patches import Circle, Rectangle
		if view == 'side':
			return Rectangle((cyl.center[0]-cyl.radius, cyl.center[2]-cyl.length/2.), cyl.radius*2, cyl.length, **kwargs)
		elif view == 'top':
			return Circle((cyl.center.x, cyl.center.y), radius=cyl.radius, **kwargs)
	
	def Finish(self):
		import matplotlib
		matplotlib.use('agg')
		import pylab
		
		fig = pylab.figure(figsize=(10, 4))
		fig.subplots_adjust(wspace=0.3, bottom=0.15)
		ax = pylab.subplot(1,2,1)

		for track in self.tracks[:50]:
			l = inner.intersection(track.pos, track.dir).first
			pylab.arrow(track.pos.x, track.pos.y, l*track.dir.x, l*track.dir.y, edgecolor='k', head_width=10, width=1e-2)
			pylab.scatter([track.pos.x], [track.pos.y])
		ax.add_artist(self.cylinder_patch(outer, 'top', facecolor="None", edgecolor='r'))
		ax.add_artist(self.cylinder_patch(inner, 'top', facecolor="None", edgecolor='r'))
		
		ax.set_aspect('equal')
		ax.set_xlabel('x [m]')
		ax.set_ylabel('y [m]')
		ax.set_ylim((-1000, 1000))
		ax.set_xlim((-1000, 1000))
		
		ax = pylab.subplot(1,2,2)
		
		for track in self.tracks:
			if abs(track.pos.y) < 20:
				l = inner.intersection(track.pos, track.dir).first
				pylab.arrow(track.pos.x, track.pos.z, l*track.dir.x, l*track.dir.z, edgecolor='k', head_width=10, width=1e-2)
				pylab.scatter([track.pos.x], [track.pos.z])
		ax.add_artist(self.cylinder_patch(outer, 'side', facecolor="None", edgecolor='r'))
		ax.add_artist(self.cylinder_patch(inner, 'side', facecolor="None", edgecolor='r'))
		
		ax.set_aspect('equal')
		ax.set_xlabel('x [m]')
		ax.set_ylabel('z [m]')
		ax.set_ylim((-1000, 1000))
		ax.set_xlim((-1000, 1000))
		
		
		pylab.savefig(args.outfile)
tray.AddModule(Harvest)


tray.Execute()

