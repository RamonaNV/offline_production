#!/usr/bin/env python
from icecube import icetray, clsim, phys_services, dataclasses, simclasses
from os.path import expandvars, join, isfile
from os import environ
import tempfile
import shutil
from argparse import ArgumentParser
from collections import defaultdict

import numpy

parser = ArgumentParser()
parser.add_argument('-g', '--gcd-file', default=expandvars('$I3_DATA/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz'))
parser.add_argument('--oversize', default=1, type=int)
parser.add_argument('--energy', default=1e5, type=float)
parser.add_argument('-o', '--output-file', default='hits')
parser.add_argument('-p', '--cascade-position', default=(0, 0, -400))

args = parser.parse_args()

DetectorParams = clsim.traysegments.common.setupDetector(
    GCDFile=args.gcd_file,
    DOMOversizeFactor=args.oversize,
    HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.nominal"),
     IceModelLocation=expandvars("$I3_SRC/ice-models/resources/models/spice_3.2.2-for_clsim")
    # IceModelLocation=expandvars("$I3_SRC/ice-models/resources/models/spice_bfr-dv1_complete")
    # IceModelLocation=expandvars("./spice_bfr_cranked_up")
)

rng = phys_services.I3GSLRandomService(0)

clsimer = clsim.traysegments.common.setupPropagators(rng, DetectorParams,
                                                     UseCUDA=False)[0]
clsimer_CUDA = clsim.traysegments.common.setupPropagators(rng, DetectorParams,
                                                     UseCUDA=True)[0]


print('---> clsim (ocl) granularity %d, bunch size %d' % (clsimer.workgroupSize,
                                                    clsimer.maxNumWorkitems))
print('---> clsim (CUDA) granularity %d, bunch size %d' % (clsimer_CUDA.workgroupSize,
                                                  clsimer_CUDA.maxNumWorkitems))

try:
    from math import gcd
except ImportError:
    from fractions import gcd


def lcm(a, b):
    return a*b/gcd(a, b)


granularity = int(lcm(clsimer_CUDA.workgroupSize, clsimer.workgroupSize))
maxBunchSize = min((clsimer.maxNumWorkitems, clsimer_CUDA.maxNumWorkitems))
maxBunchSize -= (maxBunchSize % granularity)
print('---> common granularity %d, bunch size %d' % (granularity,
                                                     maxBunchSize))

stepGenerator = clsim.I3CLSimLightSourceToStepConverterAsync()

stepGenerator.SetLightSourceParameterizationSeries(DetectorParams['ParameterizationList'])
stepGenerator.SetMediumProperties(DetectorParams['MediumProperties'])
stepGenerator.SetRandomService(rng)
stepGenerator.SetWlenBias(DetectorParams['WavelengthGenerationBias'])
stepGenerator.SetMaxBunchSize(maxBunchSize)
stepGenerator.SetBunchSizeGranularity(granularity)
stepGenerator.Initialize()

if isinstance(args.cascade_position,  str):
    if args.cascade_position.lower() == 'random':
        cascade_position = (numpy.random.uniform(-400, 400),
                            numpy.random.uniform(-400, 400),
                            numpy.random.uniform(-500, 500))
    elif args.cascade_position.lower() == 'random_full':
        cascade_position = (numpy.random.uniform(-600, 600),
                            numpy.random.uniform(-500, 500),
                            numpy.random.uniform(-500, 500))
    elif args.cascade_position.lower() == 'center':
        cascade_position = (0, 0, 0)
    else:
        raise ValueError("Invalid value for argument cascade_position!")
elif isinstance(args.cascade_position, tuple):
    if len(args.cascade_position) == 3:
        cascade_position = args.cascade_position
else:
    raise ValueError("Invalid value for argument cascade_position!")
print("---> Cascade position is:", cascade_position)

p = dataclasses.I3Particle()
p.type = p.EMinus
p.energy = args.energy
p.time = 0
p.pos = dataclasses.I3Position(*cascade_position)
p.dir = dataclasses.I3Direction(0, 0)

for i in range(1):
    stepGenerator.EnqueueLightSource(clsim.I3CLSimLightSource(p), 0)
stepGenerator.EnqueueBarrier()

photons = defaultdict(clsim.I3CLSimPhotonSeries)

i = 0
while True:
    steps, markers, particleHistories, barrierWasReset = \
        stepGenerator.GetConversionResultWithBarrierInfoAndMarkers()

    print('---> sending %d photons in bunch %d'
          % (sum((s.num for s in steps)), i))
    clsimer_CUDA.EnqueueSteps(steps, i)
    clsimer.EnqueueSteps(steps, i)
    i += 1

    result_CUDA = clsimer_CUDA.GetConversionResult()
    result_ocl = clsimer.GetConversionResult()

    photons['CUDA'].extend(result_CUDA.photons)
    photons['ocl'].extend(result_ocl.photons)

    n_CUDA = len(result_CUDA.photons)
    n_ocl = len(result_ocl.photons)
    print('---> got {CUDA: %d, ocl: %d} photons in bunch %d (CUDA/ocl=%.2f)'
          % (n_CUDA, n_ocl, result_ocl.identifier, n_CUDA/float(n_ocl)))

    if barrierWasReset:
        break

n_ocl = len(photons['ocl'])
n_CUDA = len(photons['CUDA'])
print('total unweighted {CUDA: %d, ocl: %d} photons (CUDA/ocl=%.3f)'
      % (n_CUDA, n_ocl, n_CUDA/float(n_ocl)))

if args.output_file:
    numpy.savez(args.output_file, **photons)
