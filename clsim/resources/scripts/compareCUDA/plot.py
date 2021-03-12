#!/usr/bin/env python
import numpy as np
from argparse import ArgumentParser
from scipy.stats import ks_2samp
import matplotlib
import matplotlib.pyplot as plt

fig_size = [10, 5]
matplotlib.rcParams.update({'figure.figsize': fig_size})

parser = ArgumentParser()
parser.add_argument('-i', '--input-file', default='hits.npz')
parser.add_argument('-o', '--output-file', default=None)
args = parser.parse_args()

print("Loading saved hits...")
photons = np.load(args.input_file, allow_pickle=True)
n_CUDA = len(photons['CUDA'])
n_ocl = len(photons['ocl'])
print('total unweighted {CUDA: %d, ocl: %d} photons (CUDA/ocl=%.3f)'
      % (n_CUDA, n_ocl, n_CUDA/n_ocl))

print("Summing up counts per DOM...")
cl_counts = np.zeros(shape=(86, 60), dtype=np.int)
for photon in photons['ocl']:
    cl_counts[photon.stringID - 1][photon.omID - 1] += 1

max_hits_om = np.unravel_index(cl_counts.argmax(), cl_counts.shape)
print("DOM with most hits (0 indexing):", max_hits_om)
print("Number of hits for that DOM:", cl_counts[max_hits_om])

# max_hits_om = (85, 47)

print("Creating time histogram for DOM with most hits...")
cl_times = []
for photon in photons['ocl']:
    if (photon.stringID - 1, photon.omID - 1) == max_hits_om:
        cl_times.append(photon.time)
cl_times = np.asarray(cl_times)

CUDA_times = []
for photon in photons['CUDA']:
    if (photon.stringID - 1, photon.omID - 1) == max_hits_om:
        CUDA_times.append(photon.time)
CUDA_times = np.asarray(CUDA_times)
n_CUDA, bins, _ = plt.hist(CUDA_times, bins='auto',
                          range=[CUDA_times.min(), CUDA_times.mean()*1.3],
                          histtype='step', label='CUDA')
n_ocl, _, _ = plt.hist(cl_times, bins=bins, histtype='step', label='ocl')
plt.xlabel("Arrival time at DOM {} on string {}".format(max_hits_om[1] + 1,
                                                        max_hits_om[0] + 1))
plt.ylabel("Number of hits")
plt.legend()
if isinstance(args.output_file, str):
    plt.savefig(args.output_file)
else:
    plt.savefig('{}_{}_{}.png'.format(cl_counts[max_hits_om], max_hits_om[0] +
                                      1, max_hits_om[1] + 1))
print("Saved as max_hist.png.")

print("KS test result:")
print(ks_2samp(n_CUDA, n_ocl))
print("\nInfo: The p-value should be at least 99% for sufficiently large "
      "number of hits in the brightest DOM.")
print("The total hits should match with a difference of only ~0.1%.")
