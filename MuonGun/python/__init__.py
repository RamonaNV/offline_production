from icecube.load_pybindings import load_pybindings
import icecube.icetray, icecube.dataclasses, icecube.simclasses, icecube.phys_services# be nice and pull in our dependencies
load_pybindings(__name__,__path__)

import inspect
def monkeypatch(cls):
    def monkeypatch(f):
        name = f.__name__
        if inspect.getargspec(f).args[0] == 'cls':
            f = classmethod(f)
        setattr(cls, name, f)
        return None
    return monkeypatch

@monkeypatch(ExtrudedPolygon)
def from_I3Geometry(cls, i3geo, padding=0):
    from collections import defaultdict
    import numpy
    strings = defaultdict(list)
    for omkey, omgeo in i3geo.omgeo:
        if omgeo.omtype != omgeo.IceTop:
            strings[omkey.string].append(list(omgeo.position))
    mean_xy = [numpy.mean(positions, axis=0)[0:2] for positions in strings.values()]
    zmax = max(max(p[2] for p in positions) for positions in strings.values())
    zmin = min(min(p[2] for p in positions) for positions in strings.values())
    
    positions = [icecube.dataclasses.I3Position(x,y,z) for x,y in mean_xy for z in (zmin,zmax)]
    
    return cls(positions, padding)
    
@monkeypatch(ExtrudedPolygon)
def from_file(cls, fname, padding=0):
    from icecube import icetray, dataio, dataclasses
    f = dataio.I3File(fname)
    fr = f.pop_frame(icetray.I3Frame.Geometry)
    f.close()
    return cls.from_I3Geometry(fr['I3Geometry'], padding)

def all_models():
	"""
	Get a list of models that are valid arguments to load_model()
	"""
	from os.path import expandvars, basename, splitext
	from glob import glob
	basedir=expandvars('$I3_BUILD/MuonGun/resources/tables/')
	suffix = '.single_flux.fits'
	return [p[:-len(suffix)] for p in map(basename, glob(basedir+'*'+suffix))]
	
def load_model(base):
	from os.path import exists, expandvars, join
	if not exists(base+'.single_flux.fits'):
		basedir=expandvars('$I3_BUILD/MuonGun/resources/tables')
		icecube.icetray.logging.log_info("The table path %s does not exist! Assuming you meant %s..." % (base, join(basedir, base)), unit="MuonGun")
		base = join(basedir, base)
	return BundleModel(SplineFlux(base+'.single_flux.fits', base+'.bundle_flux.fits'),
	    SplineRadialDistribution(base+'.radius.fits'),
	    SplineEnergyDistribution(base+'.single_energy.fits', base+'.bundle_energy.fits'))

from os.path import expandvars
def corsika_genprob(config='CascadeOptimized5Comp', atmosphere=12, hadronic_model='sibyll',
    basedir=expandvars('$I3_BUILD/MuonGun/resources/tables')):
	from os.path import join
	base = join(basedir, '%s_atmod%d_%s' % (config, atmosphere, hadronic_model.upper()))
	return CORSIKAGenerationProbability(Cylinder(1600, 800),
	    SplineFlux(base+'.single_flux.fits', base+'.bundle_flux.fits'),
	    SplineRadialDistribution(base+'.radius.fits'),
	    SplineEnergyDistribution(base+'.single_energy.fits', base+'.bundle_energy.fits'))
del expandvars
