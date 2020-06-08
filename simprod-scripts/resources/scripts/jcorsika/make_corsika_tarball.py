#!/usr/bin/env python

import subprocess
import os
import sys
import glob
import re
import shutil
import itertools
import copy

from corsika_dev.corsika_builder import CorsikaBuilder

import argparse
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--version")
args, extra_args = parser.parse_known_args()
if args.version:
    try:
        builder = CorsikaBuilder(args.version[:5]) # I add an arbitrary string after the version number when I make changes
    except Exception as e:
        print e
        exit(0)
else:
    builder = None


m = re.search('([0-9]+\.[0-9]+)\.[0-9]+', subprocess.check_output(['gfortran','--version']))
if m:
    gfortran_version = m.group(1)
else:
    print('failed doing gfortran --version')
    exit(1)

params = {
    'version': (args.version if bool(builder) else None),
    'platform': 'gfortran%s'%gfortran_version, # this 'platform' is just a tag for the tarball.
    'fluka_tarball': 'fluka2011.2c-linux-gfor64bitAA.tar.gz', # relative to src_dir
    'src_dir': os.path.join(os.environ['PWD'], 'source_tarballs'),
    'work_dir': os.environ['PWD'],
    'fortran':'gfortran', # used to specify for corsika and for fluka
    'oldnames': False
    }

corsika_options = { 'float': 'float_default',
                    'he_model': [],
                    'le_model': ['fluka'],
                    'geometry': ['flat'],
                    #'basic': ['slant', 'neutrino', 'muprod'],
                    #'basic': ['slant', 'neutrino', 'muprod', 'coast', 'icecube_fifo'],
                    'basic': [],
                    #'extra': ['thin', 'curved']
                    'extra': []
                    }

"""
make_corsika_tarball.py --version=75700 --he-model=sibyll --he-model=epos --he-model=qgsjet_ii \
                        --le-model=gheisha --le-model=fluka \
                        --basic slant \
                        --extra thin --extra neutrino --extra curved 

make_corsika_tarball.py --version=75700 --he-model=sibyll --he-model=epos --he-model=qgsjet_ii \
                        --le-model=gheisha --le-model=fluka \
                        --basic slant --basic muprod --basic history \
                        --extra thin --extra neutrino --extra curved

make_corsika_tarball.py --version=75700 --he-model=sibyll --he-model=epos --he-model=qgsjet_ii \
                        --le-model=gheisha --le-model=fluka \
                        --basic slant --basic muprod --basic history --basic icecube_output \
                        --extra thin --extra neutrino --extra curved

make_corsika_tarball.py --version=75700 --he-model=sibyll --he-model=epos --he-model=qgsjet_ii \
                        --le-model=gheisha --geometry=cylindrical \
                        --basic curved --basic icecube_output --basic charm --basic neutrino

make_corsika_tarball.py --version=75700 --platform dev --tar
"""

from optparse import OptionParser, OptionGroup
parser = OptionParser()
group = OptionGroup(parser, 'General parameters')
group.add_option("--version", help="Corsika version")
group.add_option("--platform", help="Tag for the tarball (%s)"%(params['platform']), default=params['platform'])
group.add_option("--fluka-tarball", help="Fluka tarball relative to src-dir (%s)"%(params['fluka_tarball']), default=params['fluka_tarball'])
group.add_option("--src-dir", help="directory where tarballs are located (%s)"%(params['src_dir']), default=params['src_dir'])
group.add_option("--work-dir", help="%s"%(params['work_dir']), default=params['work_dir'])
group.add_option("--fortran", help="%s"%(params['fortran']), default=params['fortran'])
group.add_option("--oldnames", help="Use old name convention (False)", default=False, action="store_true")
group.add_option("--tar", help="Create tarball", default=False, action="store_true")
parser.add_option_group(group)

if builder:
    he_models = ', '.join([o.name for o in builder.options.options.values() if o.group=='he_model'])
    le_models = ', '.join([o.name for o in builder.options.options.values() if o.group=='le_model'])
    geometries = ', '.join([o.name for o in builder.options.options.values() if o.group=='geometry'])
else:
    he_models = ''
    le_models = ''
    geometries = ''

group = OptionGroup(parser, 'CORSIKA options')
group.add_option("--float", help="Float option (%s)"%(corsika_options['float']),
                 default=corsika_options['float'])
group.add_option("--he-model", action='append',
                 help="High-energy model (%s). At least one of\n:  %s"%(', '.join(corsika_options['he_model']), he_models))
group.add_option("--le-model", action='append',
                 help="Low-energy model (%s). At least one of\n:  %s"%(', '.join(corsika_options['le_model']), le_models))
group.add_option("--geometry", action='append',
                 help="Geometry (%s). At least one of: %s"%(', '.join(corsika_options['geometry']), geometries))
group.add_option("--basic", action='append',
                 help="Basic options. All executables use them (%s)"%(corsika_options['basic']))
group.add_option("--extra", action='append',
                 help="Extra options. All combinations are used (%s)"%(corsika_options['extra']))
parser.add_option_group(group)

ops, args = parser.parse_args()
ops = vars(ops)

if not ops['version']:
    print 'Error: need to specify version\n'
    parser.print_help()
    parser.exit()

for k in params.keys():
    if not ops[k] is None: params[k] = ops[k]
for k in corsika_options.keys():
    if not ops[k] is None: corsika_options[k] = ops[k]

params['fluka_tarball'] = os.path.abspath(os.path.join(params['src_dir'], params['fluka_tarball']))
os.environ['F77']=params['fortran']

#print corsika_options

log_file = open("corsika_tarball.log", "w")

def build():
    # this is a hack so it skips when there is nothing to build
    if not corsika_options['he_model'] or not corsika_options['le_model']: return

    if not os.path.exists('{src_dir}/corsika-{version}.tar.gz'.format(**params)):
        print('CORSIKA tarball {src_dir}/corsika-{version}.tar.gz does not exist'.format(**params))
        exit(1)

    # build Fluka
    # check if tarball exists here
    # and also check if it builds (old compiler in cobalt, for example)
    if corsika_options['le_model'] == 'fluka' or 'fluka' in corsika_options['le_model']:
        if float(gfortran_version) < 4.6:
            print('gfortran version >= 4.6 is needed for Fluka')
            exit(1)
        print "building Fluka"
        if corsika_options['float'] == 'float_32b':
            params['fluka_tarball'].replace('gfor64bit','')
        fluka_dir = os.path.join(params['work_dir'], os.path.basename(params['fluka_tarball']).replace('.tar.gz',''))
        if not os.path.exists(fluka_dir) or not glob.glob(fluka_dir + '/libflukahp.a'):
            builder.build_fluka(fluka_dir, params['fluka_tarball'], params['fortran'], log_file=log_file)
        os.environ['FLUPRO']=os.path.abspath(fluka_dir)

    # build corsika
    os.chdir(params['work_dir'])
    #if os.path.exists('corsika-{version}'.format(**params)):
    #    shutil.rmtree('corsika-{version}'.format(**params))

    if not os.path.exists('corsika-{version}'.format(**params)):
        subprocess.call(['tar', 'xzf', '{src_dir}/corsika-{version}.tar.gz'.format(**params)], stdout=log_file, stderr=log_file)
    os.chdir('corsika-{version}'.format(**params))

    for geometry in corsika_options['geometry']:
        for le_model in corsika_options['le_model']:
            for he_model in corsika_options['he_model']:
                for i in range(len(corsika_options['extra'])+1)[::-1]: # reverse order!
                    # the reason for the reverse order is that, with the old names,
                    # you might create a new executable, removing an existing one, and then rename it
                    for extra in itertools.combinations(corsika_options['extra'], i):
                        options = copy.deepcopy(corsika_options)
                        options['geometry'] = geometry
                        options['he_model'] = he_model
                        options['le_model'] = le_model
                        options['other'] = options['basic'] + list(extra)
                        opt = builder.options.get_coconut_options(options)
                        if params['oldnames']:
                            name = builder.options.old_name_from_options(opt)
                        else:
                            name = builder.options.name_from_options(opt)
                        #print "building CORSIKA %s"%name
                        if glob.glob('run/corsika*' + name):
                            options['other'] = ' '.join(options['other'])
                            print 'building {he_model} {le_model} {float} {geometry} {other}'.format(**options)
                            continue
                        executable = builder.coconut(opt,
                                                     log_file,
                                                     suffix=name)
                        options['other'] = ' '.join(options['other'])
                        opt['other'] = ' '.join(opt['other'])
                        print 'building {he_model} {le_model} {float} {geometry} {other}'.format(**options)

    f = open('corsika_options', 'w')
    for executable in glob.glob('run/corsika*'):
        executable = os.path.basename(executable)
        name = executable[executable.find('_')+1:]
        f.write(executable + '\n')
        opt = builder.options.options_from_name(name, names=True)
        opt['other'] = ' '.join(opt['other'])
        f.write('   {he_model} {le_model} {float} {geometry} {other}\n'.format(**opt))
        #opt = builder.options.options_from_name(name, names=False)
        #opt['other'] = ' '.join(opt['other'])
        #f.write('   {he_model} {le_model} {float} {geometry} {other}\n'.format(**opt))
    f.close()

def tarball():
    # make tarball
    print "making tarball"
    os.chdir(params['work_dir'])
    if os.path.exists('corsika-v{version}.{platform}'.format(**params)):
        shutil.rmtree('corsika-v{version}.{platform}'.format(**params))
    os.makedirs('corsika-v{version}.{platform}'.format(**params))

    shutil.copy('corsika-{version}/corsika_options'.format(**params),
                'corsika-v{version}.{platform}/corsika_options'.format(**params))
    shutil.copy(os.path.join(os.path.dirname(__file__), 'corsika_env.sh'),
                'corsika-v{version}.{platform}'.format(**params))
    shutil.copytree('corsika-{version}/run'.format(**params),
                    'corsika-v{version}.{platform}/bin'.format(**params))
    if corsika_options['le_model'] == 'fluka' or 'fluka' in corsika_options['le_model']:
        fluka_dir = os.path.join(params['work_dir'], os.path.basename(params['fluka_tarball']).replace('.tar.gz',''))
        shutil.copytree(fluka_dir,
                        'corsika-v{version}.{platform}/fluka'.format(**params))
    os.makedirs('corsika-v{version}.{platform}/lib'.format(**params))

    epos_files = [f for f in glob.glob('corsika-{version}/epos/*.*'.format(**params))
                  if not re.search('.f$', f) and not re.match('Make', f) and not re.match('README', f)]
    for f in epos_files:
        shutil.copy(f, 'corsika-v{version}.{platform}/bin'.format(**params))

    #copy libg2c.so to the lib directory (only if required)
    fluka_executables = [f for f in glob.glob('corsika-{version}/run/corsika*_fluka*'.format(**params)) if os.access(f, os.X_OK)]
    if fluka_executables:
        ldd = subprocess.check_output(['ldd', fluka_executables[0]], stderr=log_file)
        if re.search('libg2c', ldd):
            libs = [l.strip().split() for l in ldd.split('\n') if re.search('libg2c', l)]
            if libs:
                shutil.copy(libs[0][2], 'corsika-v{version}.{platform}/lib'.format(**params))
            else:
                print "ldd failed:\n", ldd
                sys.exit(0)


    # copy COAST library if it is there
    if 'COAST_USER_LIB' in os.environ and ('coast' in corsika_options['basic'] or 'coast' in corsika_options['extra']):
        for l in glob.glob(os.path.join(os.environ['COAST_USER_LIB'], 'lib*.so')):
            shutil.copy(l, 'corsika-v{version}.{platform}/lib'.format(**params))

    print "finishing tarball"
    subprocess.call(['tar', 'czf', 'corsika-v{version}.{platform}.tar.gz'.format(**params), 'corsika-v{version}.{platform}'.format(**params)], stdout=log_file, stderr=log_file)
    f = open('corsika-v{version}.{platform}.md5sum'.format(**params), 'w')
    subprocess.call(['md5sum', 'corsika-v{version}.{platform}.tar.gz'.format(**params)], stdout=f)
    f.close()

    # make sure QGSJet and EPOS data is there


build()
if ops['tar']: tarball()
