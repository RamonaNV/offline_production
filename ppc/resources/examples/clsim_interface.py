#!/usr/bin/env python

"""
Demo: feeding the same steps through CLSim and PPC
"""

from __future__ import print_function

from icecube import icetray, ppc, clsim, phys_services, dataclasses, simclasses
from os.path import expandvars, join, isfile
from os import environ
import tempfile, shutil

import numpy

parser = ArgumentParser()
parser.add_argument('-g', '--gcd-file', default=expandvars('$I3_DATA/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz'))
parser.add_argument('--use-gpus', default=False, action='store_true')
parser.add_argument('--oversize', default=1, type=int)
parser.add_argument('--energy', default=1e3, type=float)
parser.add_argument('-o', '--output-file', default=None)

args = parser.parse_args()

DetectorParams = clsim.traysegments.common.setupDetector(
    GCDFile=args.gcd_file,
    DOMOversizeFactor=args.oversize,
    HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.nominal"),
    IceModelLocation=expandvars("$I3_SRC/ice-models/resources/models/spice_3.2.2-for_clsim")
)
icetray.logging.set_level_for_unit('ppc', 'WARN')

rng = phys_services.I3GSLRandomService(0)

from icecube.ppc import MakeCLSimPropagator
ppcer = MakeCLSimPropagator(DetectorParams, UseGPUs=args.use_gpus, UseCPUs=not args.use_gpus)

clsimer = clsim.traysegments.common.setupPropagators(rng, DetectorParams, UseCPUs=not args.use_gpus)[0]

print('---> ppc granularity %d, bunch size %d' % (ppcer.workgroupSize, ppcer.maxNumWorkitems))
print('---> clsim granularity %d, bunch size %d' % (clsimer.workgroupSize, clsimer.maxNumWorkitems))

try:
    from math import gcd
except ImportError:
    from fractions import gcd
lcm = lambda a,b: a*b/gcd(a,b)
granularity = int(lcm(ppcer.workgroupSize, clsimer.workgroupSize))
maxBunchSize = min((clsimer.maxNumWorkitems, ppcer.maxNumWorkitems))
maxBunchSize -= (maxBunchSize % granularity)
print('---> common granularity %d, bunch size %d' % (granularity, maxBunchSize))

stepGenerator = clsim.I3CLSimLightSourceToStepConverterAsync()

stepGenerator.SetLightSourceParameterizationSeries(DetectorParams['ParameterizationList'])
stepGenerator.SetMediumProperties(DetectorParams['MediumProperties'])
stepGenerator.SetRandomService(rng)
stepGenerator.SetWlenBias(DetectorParams['WavelengthGenerationBias'])
stepGenerator.SetMaxBunchSize(maxBunchSize)
stepGenerator.SetBunchSizeGranularity(granularity)
stepGenerator.Initialize()

p = dataclasses.I3Particle()
p.type = p.EMinus
p.energy = args.energy
p.time = 0
p.pos = dataclasses.I3Position(0,0,-400)
p.dir = dataclasses.I3Direction(0,0)

for i in range(1):
	stepGenerator.EnqueueLightSource(clsim.I3CLSimLightSource(p), 0)
stepGenerator.EnqueueBarrier()

from collections import defaultdict
photons = defaultdict(clsim.I3CLSimPhotonSeries)

i = 0
while True:
	steps, markers, particleHistories, barrierWasReset = stepGenerator.GetConversionResultWithBarrierInfoAndMarkers()

	print('---> sending %d photons in bunch %d' % (sum((s.num for s in steps)), i))
	ppcer.EnqueueSteps(steps, i)
	clsimer.EnqueueSteps(steps, i)
	i += 1
	
	result_ppc =  ppcer.GetConversionResult()
	result_clsim =  clsimer.GetConversionResult()
	
	photons['ppc'].extend(result_ppc.photons)
	photons['clsim'].extend(result_clsim.photons)

	n_ppc = len(result_ppc.photons)
	n_clsim = len(result_clsim.photons)
	print('---> got {ppc: %d, clsim: %d} photons in bunch %d (ppc/clsim=%.2f)' % (n_ppc, n_clsim, result_ppc.identifier, n_ppc/float(n_clsim)))
	
	if barrierWasReset:
		break

n_ppc = len(photons['ppc'])
n_clsim = len(photons['clsim'])
print('total unweighted {ppc: %d, clsim: %d} photons (ppc/clsim=%.3f)' % (n_ppc, n_clsim, n_ppc/float(n_clsim)))
n_ppc = sum([p.weight for p in photons['ppc']])
n_clsim = sum([p.weight for p in photons['clsim']])
print('total weighted {ppc: %d, clsim: %d} photons (ppc/clsim=%.3f)' % (n_ppc, n_clsim, n_ppc/float(n_clsim)))

# Down-convert to MCPEs. If the Cherenkov photon generator is configured
# correctly, this should remove any constant or wavelength-dependent scale
# factors that may be present in the raw number of detected photons.
hitter = clsim.I3CLSimPhotonToMCPEConverterForDOMs(
    rng,
    DetectorParams['WavelengthAcceptance'],
    DetectorParams['AngularAcceptance']
)
def generate_mcpe(photons):
	for clsim_photon in photons:
		p = simclasses.I3CompressedPhoton()
		p.dir = clsim_photon.dir
		p.pos = clsim_photon.pos
		p.time = clsim_photon.time
		p.wavelength = clsim_photon.wavelength
		p.weight = clsim_photon.weight
		key = dataclasses.ModuleKey(clsim_photon.stringID, clsim_photon.omID)
		hit = hitter.Convert(key, p)
		if hit:
			yield hit[1]

n_ppc = len(list(generate_mcpe(photons['ppc'])))
n_clsim = len(list(generate_mcpe(photons['clsim'])))
print('total {ppc: %d, clsim: %d} MCPE (ppc/clsim=%.3f)' % (n_ppc, n_clsim, n_ppc/float(n_clsim)))

if args.output_file:
    numpy.savez(args.output_file, **photons)
