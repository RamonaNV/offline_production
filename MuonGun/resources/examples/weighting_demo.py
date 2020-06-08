#!/usr/bin/env python

from argparse import ArgumentParser
from os.path import expandvars
parser = ArgumentParser()
parser.add_argument("infiles", nargs="+", help="MuonGun input files (e.g. from generate_single_muons.py)")
parser.add_argument("-o", "--output-hdf", default="muongun_weights.hdf5")
parser.add_argument("outfile", help="output plot file")
args = parser.parse_args()

from icecube import icetray, dataclasses, dataio, MuonGun, phys_services
from I3Tray import I3Tray

icetray.logging.set_level_for_unit("MuonGun", "INFO")
icetray.logging.set_level_for_unit("I3MuonGun::WeightCalculator", "TRACE")

def harvest_generators(infiles):
	"""
	Harvest serialized generator configurations from a set of I3 files.
	"""
	from icecube.icetray.i3logging import log_info as log
	generator = None
	for fname in infiles:
		f = dataio.I3File(fname)
		fr = f.pop_frame(icetray.I3Frame.Stream('S'))
		f.close()
		if fr is not None:
			for k in fr.keys():
				v = fr[k]
				if isinstance(v, MuonGun.GenerationProbability):
					log('%s: found "%s" (%s)' % (fname, k, type(v).__name__), unit="MuonGun")
					if generator is None:
						generator = v
					else:
						generator += v
	return generator

def book_weights(infiles, outfile, model='Hoerandel5_atmod12_SIBYLL'):
	tray = I3Tray()
	
	tray.AddModule('I3Reader', 'reader', filenamelist=infiles)
	
	tray.AddModule(lambda frame: frame.Put('MCPrimary', frame['I3MCTree'].primaries[0]), Streams=[icetray.I3Frame.DAQ])
	tray.AddModule(lambda frame: frame.Put('Muon', frame['I3MCTree'][1]), Streams=[icetray.I3Frame.DAQ])

	model = MuonGun.load_model(model)
	generator = harvest_generators(infiles)
	tray.AddModule('I3MuonGun::WeightCalculatorModule', 'MuonWeight', Model=model, Generator=generator)
	from icecube.hdfwriter import I3SimHDFWriter
	tray.AddSegment(I3SimHDFWriter, 'scribe',
	    Output=outfile,
	    Keys=['MCPrimary', 'Muon', 'MuonWeight',
	        dict(key='I3MCTree', name='BundleParameters',
	             converter=MuonGun.converters.MuonBundleConverter(1, generator.surface))],
	    Types=[]
	)
	
	
	tray.Execute()
	

book_weights(args.infiles, args.output_hdf)

try:
	import matplotlib
	matplotlib.use('agg')
	import tables, pylab, numpy
	
	with tables.open_file(args.output_hdf) as hdf:
		axis = hdf.root.MCPrimary.read()
		bundle = hdf.root.BundleParameters.read()	
		# We can get pre-calculated weights calcuated by the icetray module
		muon_energy = hdf.root.Muon.cols.energy[:]
		pre_weights = hdf.root.MuonWeight.cols.value[:]
	
	# or make a freestanding WeightCalculator with a different flux model
	model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
	generator = harvest_generators(args.infiles)
	weighter = MuonGun.WeightCalculator(model, generator)
	
	# and use the bundle axis and its radius/energy distribution to calculate a weight
	try:
		post_weights = weighter(axis['x'], axis['y'], axis['z'], axis['zenith'], axis['azimuth'],
		    bundle['multiplicity'], bundle['energy'], bundle['radius'])
	except:
		# work around lack of vectorized pybindings
		@numpy.vectorize
		def weight(x,y,z,zenith,azimuth,multiplicity,e,r):
			axis = dataclasses.I3Particle()
			axis.pos = dataclasses.I3Position(x,y,z)
			axis.dir = dataclasses.I3Direction(zenith,azimuth)
			assert multiplicity == 1
			bundle = MuonGun.BundleConfiguration()
			bundle.append(MuonGun.BundleEntry(float(r),float(e)))
			return weighter(axis, bundle)
		post_weights = weight(axis['x'], axis['y'], axis['z'], axis['zenith'], axis['azimuth'],
		    bundle['multiplicity'], bundle['energy'], bundle['radius'])
	
	pylab.hist(muon_energy, bins=numpy.logspace(4, 7, 101), log=True, histtype='bar', label='unweighted (counts)')
	pylab.hist(muon_energy, weights=post_weights, bins=numpy.logspace(4, 7, 101), log=True, histtype='bar', label='GaisserH4a')
	pylab.hist(muon_energy, weights=pre_weights, bins=numpy.logspace(4, 7, 101), log=True, histtype='bar', label='Hoerandel')
	
	pylab.legend(loc='lower left')
	pylab.xlabel('Muon energy at sampling surface')
	pylab.ylabel('Rate [Hz]')
	pylab.loglog()
	
	pylab.savefig(args.outfile)
except ImportError as e:
	print(e)
	pass

