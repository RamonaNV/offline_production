
#!/usr/bin/env python

"""
Demo: feeding the same steps through CLSim 
"""

from __future__ import print_function

from icecube import icetray, ppc, clsim, phys_services, dataclasses, simclasses
from os.path import expandvars, join, isfile
from os import environ
import tempfile, shutil
import argparse
import numpy

"""The MIT License (MIT)

Copyright (c) 2020, Ramona Hohl, rhohl@nvidia.com

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


parser =argparse.ArgumentParser()
parser.add_argument('-g', '--gcd-file', default=expandvars('$I3_TESTDATA/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz'))
parser.add_argument('--use-gpus', default=True, action='store_true')
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
 

clsimer = clsim.traysegments.common.setupPropagators(rng, DetectorParams, UseCPUs=not args.use_gpus)[0]

 
#pprint('---> clsim granularity %d, bunch size %d' % (clsimer.workgroupSize, clsimer.maxNumWorkitems))

try:
    from math import gcd
except ImportError:
    from fractions import gcd
lcm = lambda a,b: a*b/gcd(a,b)
granularity = int( clsimer.workgroupSize)
maxBunchSize = clsimer.maxNumWorkitems 
maxBunchSize -= (maxBunchSize % granularity)
#pprint('---> common granularity %d, bunch size %d' % (granularity, maxBunchSize))

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

	#print('---> sending %d photons in bunch %d' % (sum((s.num for s in steps)), i))
	 
	clsimer.EnqueueSteps(steps, i)
	i += 1
	
	 
	result_clsim =  clsimer.GetConversionResult()
	
 
	photons['clsim'].extend(result_clsim.photons)

	 
	n_clsim = len(result_clsim.photons)
#p	print('---> got  clsim: %d photons in bunch %d' % ( n_clsim, result_clsim.identifier))
	
	if barrierWasReset:
		break


n_clsim = len(photons['clsim'])
#pprint('total unweighted clsim: %d photons  ' %   n_clsim )

n_clsim = sum([p.weight for p in photons['clsim']])


#pprint('total weighted clsim: %d photons  ' %   n_clsim )

 
 
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

 
n_clsim = len(list(generate_mcpe(photons['clsim'])))
print('total clsim: %d' % ( n_clsim ))

if args.output_file:
    numpy.savez(args.output_file, **photons)
