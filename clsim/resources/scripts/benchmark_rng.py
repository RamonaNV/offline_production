
import itertools
rng_types = ['mwc'] + ['{}{}x{}'.format(*combo) for combo in itertools.product(('threefry', 'philox'), [2,4], [32,64])]
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("--rngtype", default="mwc", choices=['all'] + rng_types)
parser.add_argument("--gpu", default=False, action="store_true")
parser.add_argument("--verbose", default=False, action="store_true")

args = parser.parse_args()


from icecube import icetray, clsim, phys_services
from icecube.clsim.traysegments.common import configureOpenCLDevices
if args.verbose:
    icetray.logging.set_level('DEBUG')
else:
    icetray.logging.set_level('FATAL')

openCLDevice = configureOpenCLDevices(UseCPUs=not args.gpu, UseGPUs=args.gpu)[0]

repetitions = 10
openCLDevice.useNativeMath=True
workgroupSize = 256 if args.gpu else 1024
workItemsPerIteration = workgroupSize*128*(4 if args.gpu else 1)
iterations = 10000 if args.gpu else 16
rng = phys_services.I3SPRNGRandomService(seed=3244, nstreams=2, streamnum=0)

if args.rngtype != 'all':
    rng_types = [args.rngtype]

for rng_type in rng_types:

    tester = clsim.I3CLSimRandomNumberGeneratorBenchmark(device=openCLDevice,
                                                   workgroupSize=workgroupSize,
                                                   workItemsPerIteration=workItemsPerIteration,
                                                   iterations=iterations,
                                                   randomService=rng,
                                                   rngType="" if rng_type == "mwc" else rng_type)
    values, time = tester.GenerateRandomNumbers(repetitions)

    total_reps = workItemsPerIteration*iterations*repetitions
    print("{}: {:.4f} ns/rep".format(rng_type, time/float(total_reps)))
