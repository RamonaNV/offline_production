
from cubicle.dagman import DAG, Node
import glob
import os, shutil

hmodel = "SIBYLL"
atmod = 12
outdir = "/data/uwa/jvansanten/projects/2012/muongun/corsika/%(hmodel)s/atmod%(atmod)s" % locals()
os.mkdir(outdir)

filesets = dict()
basedir = "/data/sim/IceCube/2011/generated/dcorsika/length_1600_radius_800_v6960-5comp_sibyll_5-component/emin_600_emax_1e11_dslope_0_pgam_E2.0E2.0E2.0E2.0E2.0"
print('%(basedir)s/atmod_%(atmod)d/*/corsika_2500000_*.gz' % locals())
filesets['standard'] = dict(
    files=sorted(glob.glob('%(basedir)s/atmod_%(atmod)d/*/corsika_2500000_*.gz' % locals())),
    # chunk=100,
    chunk=1,
    
)
basedir = "/data/sim/IceCube/2011/generated/dcorsika/9231/length_1600_radius_800_v6960-5comp_sibyll_5-component/emin_5e4_emax_1e11_dslope_0_pgam_E2.6E2.6E2.6E2.6E2.6"
filesets['he'] = dict(
    files=sorted(glob.glob('%(basedir)s/atmod_%(atmod)d/*/corsika_1000000_*.gz' % locals())),
    # chunk=20,
    chunk=1,
)

weight_args = " ".join(['--n-%s=%d' % (k, len(fs['files'])) for k, fs in filesets.items()])
print('weight_args: ', weight_args)

def chunks(l, n):
	for i in range(0, len(l), n):
		yield l[i:i+n]

dag = DAG()
i = 0;
jobs = []
for (label, fileset) in list(filesets.items()):
	for files in chunks(fileset['files'], fileset['chunk']):
		i += 1
		paths = " ".join(files)
		job = 'corsika_%(hmodel)s_atmod%(atmod)d_%(i).6d.hdf5' % locals()
		args = "%(paths)s %(outdir)s/%(job)s %(weight_args)s" % locals() 
		node = Node(job, 'fill_histograms.sub', args=args)
		dag.add(node)
		jobs.append(node)

paths = " ".join([outdir+'/'+node.label for node in jobs])
node = Node('corsika_%(hmodel)s_atmod%(atmod)s_merge' % locals(), 'histadd.sub', arguments="--overwrite --norm=1 %(paths)s %(outdir)s.hdf5" % locals())
dag.add(node, jobs)

dag.write()
