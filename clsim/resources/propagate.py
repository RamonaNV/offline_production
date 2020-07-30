#!/usr/bin/env python
from icecube import icetray, ppc, clsim, phys_services, dataclasses, simclasses
from icecube.ppc import MakeCLSimPropagator
from os.path import expandvars, join, isfile
from os import environ
import tempfile
import shutil
from argparse import ArgumentParser
from collections import defaultdict

import numpy

parser = ArgumentParser()
parser.add_argument('-g', '--gcd-file', default=expandvars('$I3_DATA/GCD/GeoCalibDetectorStatus_IC86_Merged.i3.gz'))
parser.add_argument('--use-gpus', default=True, action='store_true')
# DOM oversizing still gives different results (PPC > clsim).
# TODO: I should track this down!
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
icetray.logging.set_level_for_unit('ppc', 'WARN')

rng = phys_services.I3GSLRandomService(0)

ppcer = MakeCLSimPropagator(DetectorParams, UseGPUs=args.use_gpus,
                            UseCPUs=not args.use_gpus) #, CopyConfig=True)

clsimer = clsim.traysegments.common.setupPropagators(rng, DetectorParams,
                                                     UseCPUs=not
                                                     args.use_gpus)[0]

print('---> ppc granularity %d, bunch size %d' % (ppcer.workgroupSize,
                                                  ppcer.maxNumWorkitems))
print('---> clsim granularity %d, bunch size %d' % (clsimer.workgroupSize,
                                                    clsimer.maxNumWorkitems))

try:
    from math import gcd
except ImportError:
    from fractions import gcd


def lcm(a, b):
    return a*b/gcd(a, b)


granularity = int(lcm(ppcer.workgroupSize, clsimer.workgroupSize))
maxBunchSize = min((clsimer.maxNumWorkitems, ppcer.maxNumWorkitems))
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
    ppcer.EnqueueSteps(steps, i)
    clsimer.EnqueueSteps(steps, i)
    i += 1

    result_ppc = ppcer.GetConversionResult()
    result_clsim = clsimer.GetConversionResult()

    photons['ppc'].extend(result_ppc.photons)
    photons['clsim'].extend(result_clsim.photons)

    n_ppc = len(result_ppc.photons)
    n_clsim = len(result_clsim.photons)
    print('---> got {ppc: %d, clsim: %d} photons in bunch %d (ppc/clsim=%.2f)'
          % (n_ppc, n_clsim, result_ppc.identifier, n_ppc/float(n_clsim)))

    if barrierWasReset:
        break

n_ppc = len(photons['ppc'])
n_clsim = len(photons['clsim'])
print('total unweighted {ppc: %d, clsim: %d} photons (ppc/clsim=%.3f)'
      % (n_ppc, n_clsim, n_ppc/float(n_clsim)))

if args.output_file:
    numpy.savez(args.output_file, **photons)