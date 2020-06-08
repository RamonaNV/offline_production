"""
Very primitive test for now.
"""

def test_corsika_stage():
    import os
    import corsika_dev.corsika_stage
    par = {'url':'http://prod-exe.icecube.wisc.edu', #'url':'file:/data/sim/sim-new/downloads',
           'version':'v73700',
           'platform':'gfortran_4.8.2_IT'}

    r = corsika_dev.corsika_stage.corsika_stage(**par)
    print 'corsika-{version}.{platform} exists:'.format(**par), os.path.exists('corsika-{version}.{platform}'.format(**par))
    print 'corsika-{version}.{platform}.tar.gz exists:'.format(**par), os.path.exists('corsika-{version}.{platform}.tar.gz'.format(**par))

# the types of the parameters are not enforced yet
parameters = [{'steering':True, 'version': '73700', 'he_model': 'sibyll', 'le_model': 'fluka', 'nshow': 100,
               'prmpar': 2652, 'emin': 1e5, 'emax': 1e5, 'corsika_path': '/data/user/jgonzalez/corsika/corsika-v73700.gfortran_4.8.2_IT',
               'float': 'float_default', 'geometry': 'flat', 'old_names':True}]

command = [['./scripts/corsika.py', '--steering', '--version=73700', '--he-model=sibyll', '--le-model=fluka', '--nshow=100',
           '--prmpar=2652', '--emin=1e5', '--emax=1e5', '--corsika-path=/data/user/jgonzalez/corsika/corsika-v73700.gfortran_4.8.2_IT',
           '--float=float_default', '--geometry=flat', '--old-names']]

steering = ["""
DIRECT /data/user/jgonzalez/IceCubeSoftware/Sandbox/corsika_dev/trunk/ 
RUNNR 1 
EVTNR 1 
NSHOW 100 
PRMPAR 2652 
ESLOPE -2 
ERANGE 100000.000000 100000.000000 
THETAP 0.000000 65.000000 
PHIP -180.000000 180.000000 
SEED 3 0 0 
SEED 4 0 0 
SEED 5 0 0 
OBSLEV 283700 
ELMFLG T T 
RADNKG 2.E5 
MAGNET  16.40 -53.40 
LONGI T 20.000000 T F 
HADFLG  0  1  0  1  0  2 
ECUTS  0.1000 0.1000 0.0003 0.0003 
MAXPRT 0 
ECTMAP 100 
STEPFC 1.0 
DEBUG F 6 F 1000000 
SIBYLL  T  0 
SIBSIG  T 
ATMOD   13                         october atmosphere
EXIT 
"""]

def test_corsika_config_cmd():
    import subprocess
    for i, (cmd, out) in enumerate(zip(command, steering)):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output = p.communicate()[0]
        print 'test_corsika_config_cmd', i, output.strip()==out.strip()

def test_corsika_config_py():
    import corsika_dev
    for i, (par, out) in enumerate(zip(parameters, steering)):
        corsika = corsika_dev.CorsikaOptions(version=par['version'], defaults=par)
        output = corsika.steering()
        print 'test_corsika_config_py', i, output.strip()==out.strip()

test_corsika_stage()

test_corsika_config_cmd()
test_corsika_config_py()

