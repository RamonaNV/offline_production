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
n_ppc = len(photons['ppc'])
n_clsim = len(photons['clsim'])
print('total unweighted {ppc: %d, clsim: %d} photons (ppc/clsim=%.3f)'
      % (n_ppc, n_clsim, n_ppc/n_clsim))

print("Summing up counts per DOM...")
cl_counts = np.zeros(shape=(86, 60), dtype=np.int)
for photon in photons['clsim']:
    cl_counts[photon.stringID - 1][photon.omID - 1] += 1

max_hits_om = np.unravel_index(cl_counts.argmax(), cl_counts.shape)
print("DOM with most hits (0 indexing):", max_hits_om)
print("Number of hits for that DOM:", cl_counts[max_hits_om])

# max_hits_om = (85, 47)

print("Creating time histogram for DOM with most hits...")
cl_times = []
for photon in photons['clsim']:
    if (photon.stringID - 1, photon.omID - 1) == max_hits_om:
        cl_times.append(photon.time)
cl_times = np.asarray(cl_times)

ppc_times = []
for photon in photons['ppc']:
    if (photon.stringID - 1, photon.omID - 1) == max_hits_om:
        ppc_times.append(photon.time)
ppc_times = np.asarray(ppc_times)
n_ppc, bins, _ = plt.hist(ppc_times, bins='auto',
                          range=[ppc_times.min(), ppc_times.mean()*1.3],
                          histtype='step', label='ppc')
n_clsim, _, _ = plt.hist(cl_times, bins=bins, histtype='step', label='clsim')
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
print(ks_2samp(n_ppc, n_clsim))
print("\nInfo: The p-value should be at least 99% for sufficiently large "
      "number of hits in the brightest DOM.")
print("The total hits should match with a difference of only ~0.1%.")