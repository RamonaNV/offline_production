#!/bin/env python
#
"""
 Interface for configuring pre/post icetray scripts

 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Juan Carlos Diaz Velez <juancarlos@icecube.wisc.edu>
"""

import os
import re
import sys
import math
import glob
import time
import string
import shutil
from .. import ipmodule
from os import system
from os.path import expandvars
import logging
import getpass
try:
    from commands import getstatusoutput
except ImportError:
    from subprocess import getstatusoutput

import base64

if sys.version_info[0] >= 3:
    import urllib.parse as urlparse
    from urllib.error import HTTPError
    from urllib.request import urlopen, Request
else:
    import urlparse
    from urllib2 import HTTPError
    from urllib2 import urlopen, Request

def strip_auth(url,auth=None):
    """
    urlopen() doesn't support inline auth. Strip it out and
    construct the appropriate header by hand.
    """
    parsed = urlparse.urlparse(url)
    if '@' in parsed.netloc:
        auth, netloc = parsed.netloc.split('@')
        parts = list(parsed)
        parts[1] = netloc
        url = urlparse.ParseResult(*parts).geturl()
    req = Request(url)
    if auth is not None:
        auth = base64.encodestring(auth.encode('utf-8'))[:-1]
        req.add_header('Authorization', 'Basic ' + auth.decode('utf-8'))
    return req

def wget(url, output_path, blocksize=2**16):
    auth = None

    if url.startswith("gsiftp://"):
        dest = os.path.expandvars(output_path)
        if os.path.isdir(dest.replace('file:','')):
           dest = os.path.join(dest, os.path.basename(url))
           if not dest.startswith("file:"):
              dest = "file:" + dest
        cmd = 'globus-url-copy -cd -r -nodcau -rst-retries 5 -rst-interval 60 %s %s'% (url,dest)
        retval = os.system(cmd)
        return output_path

    for wgetrc in ('$WGETRC','$HOME/.wgetrc'):
        wgetrc = os.path.expandvars(os.path.expanduser(wgetrc))
        if os.path.exists(wgetrc):
            user = None
            passwd = None
            for line in open(wgetrc):
                parts = [x.strip() for x in line.split('=',1)]
                if parts[0] == 'http_user':
                    user = parts[1]
                elif parts[0] == 'http_passwd':
                    passwd = parts[1]
                elif parts[0] == 'http_proxy':
                    os.environ['http_proxy'] = parts[1]
            if user and passwd:
                auth = user+':'+passwd
    
    f = None
    output_file = None
    try:
        output_file = open(output_path, "wb")
        f = urlopen(strip_auth(url,auth=auth))
            
        while True:
            block = f.read(blocksize)
            output_file.write(block)
            if len(block) < blocksize:
                break
    except Exception as e:
        if os.path.exists(output_path):
            os.remove(output_path)
        raise
    finally:
        if f is not None:
            f.close()
        if output_file is not None:
            output_file.close()
    return output_path

def checksum(sum1,sum2):
   if not os.path.exists(sum1): return False
   if not os.path.exists(sum2): return False
   
   if False:
      lines1 = [x for x in open(sum1)]
      lines2 = [x for x in open(sum2)]
      diff = False
      if len(lines1) != len(lines2):
         return False
      for i in range(len(lines1)):
         m1,f1 = lines1[i].split()
         m2,f2 = lines2[i].split()
         if m1 != m2 or os.path.basename(f1) != os.path.basename(f2):
            return False
      return True
   
   if os.system("diff %s %s" % (sum1,sum2)):
      logging.debug("sums differ")
      return False
   else:
      logging.debug("sums match")
      return True

_tmpstorage = expandvars('$TMPDIR/etc/%s_icesoft' % getpass.getuser())
def fetch_tarball(meta,tmpstorage=_tmpstorage,basedir='.'):
        if ':' in tmpstorage:
            # more than one path specified, choose first good one
            storage = None
            for t in tmpstorage.split(':'):
                base = dirname(t)
                if not exists(base):
                    continue
                if not exists(t):
                    try: os.makedirs(t)
                    except: continue
                storage = t
                break
            if storage:
                tmpstorage = storage
            else:
                tmpstorage = _tmpstorage

        lockfile = os.path.join(tmpstorage,"%s.lock" % meta['filebase'])
        cache = os.path.join(tmpstorage,meta['filebase'] + meta['suffix'])
        cachedir = os.path.join(tmpstorage,meta['filebase'])
        cwd = os.getcwd()


        #fetch checksum
        md5url = meta['url'] + '.md5sum'
        logging.debug("retrieving %s..." % md5url)
        wget(md5url,os.path.join(cwd,os.path.basename(md5url)))
        # See if we need to download new tarball
        trials = 5 # limit number of attempts
        while  (not os.path.exists(cache)) or \
               (not os.path.exists(cachedir)) or \
               (not checksum(os.path.join(cwd,meta['md5sum']),os.path.join(cachedir,meta['md5sum']))): 

            os.chdir(tmpstorage)
            try: # create lockfile
                os.mknod(lockfile)
            except OSError as oserr: # file exists
                logging.error("%s %s ." % (oserr,lockfile))
                time.sleep(300) # 5min
                try:
                    if (time.time() - os.stat(lockfile).st_mtime) > 2400: 
                        os.remove(lockfile)
                except: pass
                continue # check cache again

            if os.system("rm -rf %s %s" % (cache,cachedir)):
                os.system('/usr/bin/printenv')
                raise Exception('failed to remove cache %s'%cache)
            logging.debug("retrieving %s..." % (meta['url'] + meta['suffix']))
            failed = wget(meta['url']+meta['suffix'],os.path.join(cwd,cache))
           
            if not os.path.isfile(failed):
                os.remove(lockfile)
                raise Exception("unable to fetch file from %s" % meta['url'] + meta['suffix'])

            if os.path.isdir(meta['filebase']): # clean old i3_work
                if os.system("rm -rf " + meta['metafilebase']):
                    raise Exception('failed to remove %s'%(meta['filebase']))

            # create md5sum
            os.mkdir(meta['filebase']) # create i3_work
            if os.system("tar -C%s -xzf %s" % (basedir,cache)):
                raise Exception("failed to untar %s"%cache)
            if os.system("md5sum %s > %s" % (meta['filebase']+meta['suffix'],
                                             os.path.join(cachedir,meta['md5sum']))):
                raise Exception("failed to get md5sum of %s"%(meta['filebase']+meta['suffix']))
            os.remove(lockfile)
            os.chdir(cwd)

            trials -= 1 # decrement counter
            if trials < 1: 
                raise Exception("Max no. of trails reached to stage metaproject '%s'" % meta['filebase'])

        os.chdir(cwd) # if not there already
        meta['path'] = os.path.join(tmpstorage,meta['filebase'])
        return meta

class Corsika(ipmodule.ParsingModule):
    """
    This class provides an interface for preprocessing files in iceprod
    """

    def __init__(self):
        ipmodule.ParsingModule.__init__(self)
        self.logger = logging.getLogger('iceprod::Corsika')
        self._staged = False
        self.name   = 'corsika'
        self.arrang = -119.

        self.AddParameter('version','Corsika version','v6900')
        self.AddParameter('platform','compliler platform','')
        self.AddParameter('cache','Should cache taball?',False)
        self.AddParameter('cachedir','Cache directory',
                          '$PWD/%s_icesoft/%s'% (getpass.getuser(), self.name))
        self.AddParameter('URL','fetch tarball from URL',None)
        self.AddParameter('runnum','Run number','')
        self.AddParameter('firstevent','Number of first event',1)
        self.AddParameter('seed','Random seed','1')
        self.AddParameter('egs_seed_offset','value to be added to EGS seed (for debugging)','0')
        self.AddParameter('nevents','Number of Events',0)
        self.AddParameter('outfile','Output file',"DAT%(runnum)06d") 
        self.AddParameter('logfile','Log file','%(outfile)s.log')
        self.AddParameter('outdir','Output directory','$PWD')
        self.AddParameter('topdir','Top directory','.')
        self.AddParameter('tmpdir','Temporary directory','%(topdir)s/dcors%(runnum)d')
        self.AddParameter('cvmfs', 'Path to CVMFS repository', '') 
        self.AddParameter('model','Physics Model','SIBYLL')
        self.AddParameter('lemodel','Low Energy Physics Model','gheisha')
        self.AddParameter('donkg','Run NKG',0)  # no NKG by default
        self.AddParameter('doegs','Run EGS',0)  # no EGS by default
        self.AddParameter('eslope','CR spectral index (only if ranpri=0)',-2.7)  
        self.AddParameter('crtype','CR primary type',14) 
        self.AddParameter('cthmin','Min theta of injected cosmic rays',0.0)  
        self.AddParameter('cthmax','Max theta of injected cosmic rays',89.99)  
        self.AddParameter('cphmin','Min phi of injected cosmic rays',0.0)  
        self.AddParameter('cphmax','Max phi of injected cosmic rays',360.0)  
        self.AddParameter('emin','CR min energy',600.)  
        self.AddParameter('emax','CR max energy',1.e11) 
        self.AddParameter('ecuts','hadron/em energy cut (deprecated: instead use ecuts(i),i=1..4)',0)  
        self.AddParameter('ecuts1','hadron min energy (see corsika docs)',273)  
        self.AddParameter('ecuts2','muon min energy (see corsika docs)',273)  
        self.AddParameter('ecuts3','electron min energy (see corsika docs)',0.003)  
        self.AddParameter('ecuts4','photon min energy (see corsika docs)',0.003)  
        self.AddParameter('atmod','Atmosphere model (October=13)',13) 
        self.AddParameter('nuaddi','additional information for neutrinos',False)
        self.AddParameter('obslev','distance above sea level (in cm)','2834.E2')
        self.AddParameter('flat_detector', 'Use a flat detector (ignores length and radius parameters)', False)
        self.AddParameter('length','length of generation cylinder in m (for detcfg = length/2*radius calculation)',1400.)  
        self.AddParameter('radius','radius of generation cylinder in m (for detcfg = length/2*radius calculation)',700.) 
        self.AddParameter('debug','boolean: enable disable debug mode',False) 
        self.AddParameter('dryrun','boolean: only generate INPUTS file and exit',False) 
        self.AddParameter('kcut','minimum neutrino energy required to keep the shower',None)
        self.AddParameter('ratmo','Integer: select month from 1 to 12 (Avg. profile 2007-11). 13=Avg. April 2011. Any other value would default to using the atmod parameter.',-1)
        self.AddParameter('save_long', 'Save the longitudinal profile in output file (LONG blocks).', False)
        self.AddParameter('compress', 'Compress the output file.', False)
        self.AddParameter('skipoptions', 'Compress the output file.', [])
        self.AddParameter('arrang', 'Rotation of detector relative to magnetic field direction', self.arrang) 
        self.AddParameter('curved','Using curved option?', False)

    def stage(self):
        """
        Stage files and executables
        """        
        par = self.parameters
        par['tmpdir'] = self.GetParameter('tmpdir') % par

        self.logger.info('setting up working directory: %(tmpdir)s' % par)
        if os.path.exists(par['tmpdir']):
           os.system('rm -rf %(tmpdir)s' % par)     # clean temporary directory
        if not os.path.exists(par['tmpdir']):
           os.makedirs(par['tmpdir'])     # create temporary directory
        os.chdir(par['tmpdir'])        # cd to temporary directory

        self.logger.info('caching: %s' % self.GetParameter('cache'))
        if self.GetParameter('cache'):
           cachedir = expandvars(self.GetParameter('cachedir'))
           if not os.path.exists(cachedir):
              os.makedirs(cachedir)
        else:
           cachedir = expandvars("$PWD")
        baseurl = par['url']

        meta = {}
        meta['version']  = self.GetParameter('version')
        meta['platform'] = self.GetParameter('platform')
        meta['name']     = self.name
        meta['suffix']   = ".tar.gz"
        if meta['platform']:
              meta['filebase'] = "%s-%s.%s" % (meta['name'], meta['version'],
                                               meta['platform'])
        else:
              meta['filebase'] = "%s-%s" % (meta['name'], meta['version'])
        meta['md5sum']  = "%s.md5sum" % meta['filebase']
        meta['url'] = "%s/%s" % (baseurl,meta['filebase'])

        if self.cvmfs:
           cvmfs_path = os.path.join(self.cvmfs,meta['filebase'])
           self.logger.info("CVMFS configured to '%s'",self.cvmfs)

        if self.cvmfs and os.path.exists(cvmfs_path):
           self.logger.info("CVMFS directory found: '%s'",cvmfs_path)
           par['topdir'] = cvmfs_path
        else:
           self.logger.warn("No CVMFS repo. Downloading '%s' tarball from: '%s'" %(meta['filebase'],self.url))
           fetch_tarball(meta,cachedir)
           par['topdir'] = meta['path']

        # link data files for corsika
        # Files necessary for QGSJET and QGSJET-II included
        # DPMJET and VENUS files are *not* 
        os.symlink("%(topdir)s/bin/NUCNUCCS"     % par, "NUCNUCCS")    
        os.symlink("%(topdir)s/bin/QGSDAT01"     % par, "QGSDAT01")
        os.symlink("%(topdir)s/bin/SECTNU"       % par, "SECTNU")
        os.symlink("%(topdir)s/bin/qgsdat-II-03" % par, "qgsdat-II-03")
        os.symlink("%(topdir)s/bin/sectnu-II-03" % par, "sectnu-II-03")
        os.symlink("%(topdir)s/bin/qgsdat-II-04" % par, "qgsdat-II-04")
        os.symlink("%(topdir)s/bin/sectnu-II-04" % par, "sectnu-II-04")
        os.symlink("%(topdir)s/bin/GLAUBTAR.DAT" % par, "GLAUBTAR.DAT")
        os.symlink("%(topdir)s/bin/NUCLEAR.BIN"  % par, "NUCLEAR.BIN")
        os.environ['LD_LIBRARY_PATH'] = expandvars("%(topdir)s/lib:$LD_LIBRARY_PATH"  % par)
        
        # If we have EGS files too, link those
        egsfiles = glob.glob("%(topdir)s/bin/EGSDAT*"  % par)
       
        # If we have EPOS files too, link those
        eposfiles = glob.glob("%(topdir)s/bin/epos.*"  % par)

        for file in egsfiles + eposfiles:
            os.symlink(file, os.path.basename(file))
 
        # If FLUKA exists, use it
        if os.path.exists("%(topdir)s/fluka" % par):
            os.environ['FLUPRO'] = expandvars("%(topdir)s/fluka" % par)

        #os.chdir(par['topdir'])        # cd to top directory
        return par['topdir']

        

    def Execute(self,stats):
        """ 
         Run CORSIKA
         corsika output is stdout: must create a temporary directory and cd there
        """ 
        try:
            from iceprod.core.functions import tail
        except ImportError:
            logging.error("No iceprod.core found. Will not be able to tail a file")
            def tail(*args,**kwargs):
                return ''

        cwd = os.getcwd()
        par = self.parameters

        # Retrieve tarball and stage environment
        if not self._staged:
           tmpdir = self.stage()
        self.configure()
        self.write_steering()

        # CORSIKA binary
        # New standard binary name style
        par['model'] = par['model'].upper()
        par['versionnumber'] = par['version'].lstrip('v')
        par['corsbinary'] = "%(topdir)s/bin/corsika%(versionnumber)sLinux_%(model)s_%(lemodel)s" % par
        if par['curved']:
            par['corsbinary'] += '_curved'
        if 'thinrat' in par:
            par['corsbinary'] += '_thin'
        if not os.path.exists(par['corsbinary']):
            os.chdir(cwd)           # cd back to original directory
            self.logger.error("CORSIKA binary does not exist: %(corsbinary)s\n" % par)
            self.logger.error("CORSIKA binary does not exist: corsika%(versionnumber)sLinux_%(model)s_%(lemodel)s\n" % par)
            raise Exception("CORSIKA binary does not exist: corsika%(versionnumber)sLinux_%(model)s_%(lemodel)s\n" % par)
        # Old symlink style
        #par['corsbinary'] = "%(topdir)s/bin/corsika.%(model)s.Linux" % par
        #system("cp %(corsbinary)s corsika.%(runnum)s.Linux" % par)
        os.symlink("%(corsbinary)s" % par, "corsika.%(runnum)s.Linux" % par)

        # Corsika output file 
        par['corout']  = "DAT%(runnum)06d" % par      # Plain CORSIKA output
        if par['compress']:
           par['corout']  += ".gz"

        par['logfile'] = par['logfile'] % par

        corout  = "%(outdir)s/%(corout)s" % par

        # Execution command
        cors_cmd = "%(tmpdir)s/corsika.%(runnum)d.Linux < %(inputfile)s > %(logfile)s " % par
        self.logger.info(cors_cmd)

        # Run CORSIKA
        status, output = getstatusoutput(cors_cmd)
        self.logger.info(output)
        if status: 
            os.chdir(cwd)           # cd back to original directory
            raise Exception("dCorsika python sucks! %s\n" % output)

        # check if the output is OK
        nevcorsika = 0
        try:
           status,nevcorsika=getstatusoutput("cat %(logfile)s|grep \"GENERATED EVENTS\"|grep -v NEV|awk '{ print $6 }'" % par)
           if nevcorsika.strip() != '':
               nevcorsika = int(nevcorsika.strip())
               stats['nevcorsika'] = nevcorsika
           else:
               self.logger.error('The string "GENERATED EVENTS" was not found in log file %(logfile)s' % par )
               try:
                  logfile = open(par['logfile'],'r')
                  self.logger.error(par['logfile']+':\n'+logfile.read())
                  logfile.close()
               except: pass
               status,fluka_failed=getstatusoutput('tail -n 1 %(logfile)s |grep "FLUKA TREATS LOW ENERGY HADRONIC INTERACTIONS"'%par)
               if fluka_failed:
                   self.logger.error('FLUKA initialization failed.')
        except Exception as e: 
            self.logger.error(e)
        if nevcorsika == self.GetParameter('nevents'):
            system('touch %(outdir)s/corsika.%(runnum)d.isok' % par) 
            self.logger.debug("OK")
            self.logger.debug("Corsika OK\n")
        else :
            system('touch %(outdir)s/corsika.%(runnum)d.isnotok' % par) 
            sys.stderr.write("Corsika not OK\n")
            self.logger.error("NOT OK")
            self.logger.error(tail(self.inputfilename,chars=5000))
            self.logger.error(tail(self.GetParameter('logfile'),chars=900))
            os.chdir(cwd)       # cd back to original directory
            return 1

        os.chdir(cwd)
        return 0

    def configure(self):
        """ 
         Configure and write INPUTS steering file
        """ 
        par   = self.parameters
        seed  = self.GetParameter('seed')
        egs_seed_offset = self.GetParameter('egs_seed_offset')
        par['seed1'] = int(seed)+0
        par['seed2'] = int(seed)+1+int(egs_seed_offset)
        par['seed3'] = int(seed)+2

        # NKG/EGS
        NKGparams = ""
        nkg = 'F'
        egs = 'F'
        if par['donkg']: nkg = 'T'
        if par['doegs']: egs = 'T'
        NKGparams  =    "ELMFLG  %s  %s " % (nkg,egs) 
        NKGparams +=                    "                       em. interaction flags (NKG,EGS)\n"
        if par['donkg']: # if NKG parameterizations ON
           NKGparams += "RADNKG  2.E5                           outer radius for NKG lat.dens.determ.\n"
        par['NKGparams'] = NKGparams 

        # Construct HE interaction model steering commands
        modelStr = ""
        model = self.GetParameter('model')
        if model in ("qgsjet","qgsii"):
          modelStr += "QGSJET  T  0                           use qgsjet for high energy hadrons\n"
          modelStr += "QGSSIG  T                              use qgsjet hadronic cross sections"
        elif model == "dpmjet":
          modelStr += "DPMJET  T  0                           use dpmjet for high energy hadrons\n"
          modelStr += "DPJSIG  T                              all hail Glaubtar!"
        elif model == "sibyll":
          modelStr += "SIBYLL  T  0                           use sibyll for high energy hadrons\n"
          modelStr += "SIBSIG  T                              use sibyll hadronic cross sections"
        elif model == "epos":
          modelStr += "EPOS    T  0                           use epos for high energy hadrons\n"
          modelStr += "EPOSIG  T                              use epos hadronic cross sections\n"
          modelStr += "EPOPAR input epos.param                !initialization input file for epos\n"
          modelStr += "EPOPAR fname inics epos.inics          !initialization input file for epos\n"
          modelStr += "EPOPAR fname iniev epos.iniev          !initialization input file for epos\n"
          modelStr += "EPOPAR fname initl epos.initl          !initialization input file for epos\n"
          modelStr += "EPOPAR fname inirj epos.inirj          !initialization input file for epos\n"
          modelStr += "EPOPAR fname inihy epos.ini1b          !initialization input file for epos\n"
          modelStr += "EPOPAR fname check none                !dummy output file for epos\n"
          modelStr += "EPOPAR fname histo none                !dummy output file for epos\n"
          modelStr += "EPOPAR fname data  none                !dummy output file for epos\n"
          modelStr += "EPOPAR fname copy  none                !dummy output file for epos"

        # Turn on/off dCORSIKA debugging
        if self.GetParameter('debug'):
           par['debug'] = "T"
        else:
           par['debug'] = "F"

        # Check if old-style ecuts parameter is set
        ecuts = self.GetParameter('ecuts')
        if ecuts:
           self.SetParameter('ecuts1',ecuts)
           self.SetParameter('ecuts2',ecuts)
        
        # Convert input phi from IceCube Coordinates to CORSIKA Coordinates
        # CORSIKA will rotate the particles back to IceCube Coordinates in
        # the output routine.  Also, IceCube Zenith and Azimuth point to
        # where the particle came from, while CORSIKA points to where it
        # is going.  Also CORSIKA measures zenith from -z.
        par['cphmin_cc'] = par['cphmin'] + par['arrang'] + 180
        par['cphmax_cc'] = par['cphmax'] + par['arrang'] + 180
        
        # Check the domain of min and max phi and fix if needed
        if par['cphmax_cc'] - par['cphmin_cc'] > 360.0:
          self.logger.error('Phi range greater than 360deg')
        while par['cphmin_cc'] < -360.0:
          par['cphmin_cc'] += 360.0
          par['cphmax_cc'] += 360.0
        while par['cphmax_cc'] > 360.0:
          par['cphmin_cc'] -= 360.0
          par['cphmax_cc'] -= 360.0
        
        par['nuaddi'] = self.GetParameter('nuaddi')

        length = self.GetParameter('length')
        radius = self.GetParameter('radius')
        if not self.GetParameter('flat_detector'):
            par['detcfg'] = length/(2.*radius)
        

        # Check if real atmo. is wanted
        ratmo = self.GetParameter('ratmo')

        # Write out the Corsika INPUTS steering file
        input = ""
        input += "RUNNR   %(runnum)d                        number of run\n"
        input += "EVTNR   %(firstevent)d                    number of first shower event\n"
        input += "NSHOW   %(nevents)d                       number of showers to generate\n"
        if par['crtype']:
           input += "PRMPAR  %(crtype)s                        particle type of prim. particle\n"
        input += "ESLOPE  %(eslope)f                        slope of primary energy spectrum\n"
        input += "ERANGE  %(emin)f  %(emax)f                energy range of primary particle\n"
        input += "THETAP  %(cthmin)f  %(cthmax)f            range of zenith angle (degree)\n"
        input += "PHIP    %(cphmin_cc)f  %(cphmax_cc)f            range of azimuth angle (degree)\n"
        input += "SEED    %(seed1)d   0   0                 seed for 1. random number sequence\n"
        input += "SEED    %(seed2)d   0   0                 seed for 2. random number sequence\n"
        input += "SEED    %(seed3)d   0   0                 seed for 3. random number sequence\n"
        if par['obslev']:
            input += "OBSLEV  %(obslev)s                           observation level (in cm)\n"
        input += "%(NKGparams)s"
        if par['arrang']:
            input += "ARRANG  %(arrang)f                    rotation of array to north\n"
        input += "FIXHEI  0.  0                             first interaction height & target\n"
        input += "FIXCHI  0.                                starting altitude (g/cm**2)\n"
        input += "MAGNET  16.4  -53.4                       magnetic field south pole\n"
        input += "HADFLG  0  1  0  1  0  2                  flags hadr.interact. & fragmentation\n"
        input += "%s \n" % modelStr
        input += "ECUTS   %(ecuts1).04f %(ecuts2).04f %(ecuts3).04f %(ecuts4).04f              energy cuts for particles\n"
        input += "MUADDI  T                                 additional info for muons\n"
        input += "MUMULT  T                                 muon multiple scattering angle\n"
        if par['save_long']:
            input += "LONGI   T  20.  T  F                      longit.distr. & step size & fit\n"
        else:
            input += "LONGI   F  20.  F  F                      longit.distr. & step size & fit\n"
        input += "MAXPRT  0                                 max. number of printed events\n"
        input += "ECTMAP  100                               cut on gamma factor for printout\n"
        input += "STEPFC  1.0                               mult. scattering step length fact.\n"
        input += "DEBUG   %(debug)s  6  F  1000000          debug flag and log.unit for out\n"
        #input += "DIRECT  ./                                output directory\n"
        if par['usepipe']:
           input += "PIPE T                                write to pipe\n"
        if 'detcfg' in par.keys():
            input += "DETCFG  %(detcfg)s                        detector information (l/d)\n"
        if par['nuaddi']:
            input += "NUADDI  T                             additional info for neutrinos\n"
        if par['kcut']:
            input += "KCUT    %(kcut)s                      minimum neutrino energy\n"
        #Atmosphere models created by Sam DeRidder (Gent University)
        if ratmo == 1:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -91.6956        7.01491 0.505452        -0.00181302     0.00207722\n"
           input +="ATMB        1125.71 1149.81 1032.68 490.789\n"
           input +="ATMC        821621  635444  682968  807327  5.4303203E9\n"
           input +="ATMLAY      780000  1640000 4040000 10000000\n"
        elif ratmo == 2:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -72.1988        22.7002 0.430171        -0.001203       0.00207722\n"
           input +="ATMB        1108.19 1159.77 1079.25 523.956\n"
           input +="ATMC        786271  599986  667432  780919  5.4303203E9\n"
           input +="ATMLAY      800000  1060000 4040000 10000000\n"
        elif ratmo == 3:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -63.729 -1.02799        0.324414        -0.000490772    0.00207722\n"
           input +="ATMB        1102.66 1093.56 1198.93 589.827\n"
           input +="ATMC        764831  660389  636118  734909  5.4303203E9\n"
           input +="ATMLAY      670000  2240000 4040000 10000000\n"
        elif ratmo == 4:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -69.7259        -2.79781        0.262692        -8.41695e-05    0.00207722\n"
           input +="ATMB        1111.7  1128.64 1413.98 587.688\n"
           input +="ATMC        766099  641716  588082  693300  5.4303203E9\n"
           input +="ATMLAY      760000  2200000 4040000 10000000\n"
        elif ratmo == 5:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -78.5551        -5.33239        0.312889        -9.20472e-05    0.00152236\n"
           input +="ATMB        1118.46 1169.09 1577.71 452.177\n"
           input +="ATMC        776648  626683  553087  696835  7.4095699E9\n"
           input +="ATMLAY      840000  2000000 3970000 10000000\n"
        elif ratmo == 6:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -92.6125        -8.5645 0.363986        1.65164e-05     0.00207722\n"
           input +="ATMB        1129.88 1191.98 1619.82 411.586\n"
           input +="ATMC        791177  618840  535235  692253  5.4303203E9\n"
           input +="ATMLAY      850000  1790000 3840000 10000000\n"
        elif ratmo == 7:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA	-89.9639	-13.9697	0.441631	-1.46525e-05	0.00207722\n"
           input +="ATMB	1125.73	1180.47	1581.43	373.796\n"
           input +="ATMC	784553	628042	531652	703417	5.4303203E9\n"
           input +="ATMLAY	850000	1590000	3750000	10000000\n"
        elif ratmo == 8:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -90.4253        -18.7154        0.51393 -0.00021565     0.00152236\n"
           input +="ATMB        1125.01 1175.6  1518.03 299.006\n"
           input +="ATMC        781628  633793  533269  737794  7.4095699E9\n"
           input +="ATMLAY      850000  1440000 3750000 10000000\n"
        elif ratmo == 9:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -91.686 -23.3519        0.891302        -0.000765666    0.00207722\n"
           input +="ATMB        1125.53 1169.77 1431.26 247.03\n"
           input +="ATMC        786017  645241  545022  805419  5.4303203E9\n"
           input +="ATMLAY      850000  1300000 3620000 10000000\n"
        elif ratmo == 10:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        451.616 -85.5456        2.06082 -0.001076       0.00207722\n"
           input +="ATMB        849.239 1113.16 1322.28 372.242\n"  
           input +="ATMC        225286  789340  566132  796434  5.4303203E9\n"
           input +="ATMLAY      310000  1010000 3150000 10000000\n"
        elif ratmo == 11:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -152.853        4.22741 1.38352 -0.00115014     0.00207722\n"
           input +="ATMB        1174.09 1272.49 975.906 481.615\n"
           input +="ATMC        891602  582119  643130  783786  5.4303203E9\n"
           input +="ATMLAY      850000  2240000 3240000 10000000\n"
        elif ratmo == 12:
           input += "ATMOD   10                              real atmosphere\n"
           input +="ATMA        -100.386        5.43849 0.399465        -0.00175472     0.00207722\n"
           input +="ATMB        1128.71 1198.1  858.522 480.142\n"
           input +="ATMC        829352  612649  706104  806875  5.4303203E9\n"
           input +="ATMLAY      850000  2200000 4040000 10000000\n"
        elif ratmo == 13:
           input += "ATMOD   10                              real atmosphere\n"
           input += "ATMA	-69.7259        -2.79781        0.262692        -8.41695e-05    0.00207722\n"
           input += "ATMB	1111.7  1128.64 1413.98 587.688\n"
           input += "ATMC	766099  641716  588082  693300  5.4303203E9\n"
           input += "ATMLAY	760000  2200000 4040000 10000000\n"
        else: 
            input += "ATMOD   %(atmod)s                         october atmosphere\n"

        input += "DIRECT %(outdir)s                              write ouput to given directory\n" 

        if (par['version']).strip('vV') >= "6900" and 'compress' not in par['skipoptions']:
          if par['compress']:
              input += "COMPRESS  T                          set compression setting\n" % par
          else:
              input += "COMPRESS  F                          set compression setting\n" % par

        self.steering = input


    def write_steering(self):
        # terminate steering and write it to file
        par   = self.parameters
        self.steering += "EXIT                                      terminates input\n"

        tmpdir    = par['tmpdir'] % par
        if not os.path.exists(tmpdir): os.makedirs(tmpdir)
        inputfile = open(tmpdir+"/INPUTS",'w')
        try:
           inputfile.write(self.steering % par)
        except Exception as e:
           logging.debug(self.steering)
           logging.debug(e)
           for key,opt in par.items():
               logging.debug(key, opt)
        par['inputfile'] = inputfile.name
        self.inputfilename = inputfile.name
        inputfile.close()


class dCorsika(Corsika):
    """
    This class provides an interface for preprocessing files in iceprod
    """

    def __init__(self):
        Corsika.__init__(self)
        self.name   = 'dcorsika'
        self.logger = logging.getLogger('iceprod::dCorsika')
        self.arrang = 0.

        # dcorsika specific
        self.AddParameter('ranpri','CR spectrum: 0=individual nuclei, 1=Wiebel-Sooth, 2=Hoerandel, 3= 5-component',2) 
        self.AddParameter('dslope','CR spectral index modification (only if ranpri=1,2)',0.)  
        self.AddParameter('depth','depth of the center of IceCube detector in m (for AMANDA it is 1730.)',1950.) 
        self.AddParameter('spric','separate primary energy cutoffs',True) 
        self.AddParameter('locut','Enables skew angle cutoff','T 1.58')
        self.AddParameter('pnormh','proton 5-component relative contribution',1.0)
        self.AddParameter('pnormhe','Helium 5-component relative contribution',0.1)
        self.AddParameter('pnormn','Nitrogen 5-component relative contribution',2e-3)
        self.AddParameter('pnormal','Aluminium 5-component relative contribution',6e-4)
        self.AddParameter('pnormfe','Iron 5-component relative contribution',1e-3)
        self.AddParameter('pgamh','proton 5-component spectral index',2)
        self.AddParameter('pgamhe','Helium 5-component spectral index',2)
        self.AddParameter('pgamn','Nitrogen 5-component spectral index',2)
        self.AddParameter('pgamal','Aluminium 5-component spectral index',2)
        self.AddParameter('pgamfe','Iron 5-component spectral index',2)
        self.AddParameter('f2k','Write in F2K format','T')
        
    def configure(self):
        Corsika.configure(self)

        par = self.parameters
        # dCorsika specific
        inputd  = "F2000   %(f2k)s                           choses F2000 format\n"
        inputd += "LOCUT   %(locut)s                         enables skew angle cutoff\n"
        inputd += "RANPRI  %(ranpri)s                        random primary\n"
        inputd += "SPRIC   %(spric)s                          separate primary energy cutoffs\n"
        inputd += "FSEED   F                                 enable random generator seed recovery\n"
        inputd += "DSLOPE  %(dslope)s                        slope correction\n"
        inputd += "SCURV   T 6.4E8 1.95E5                    curved surf., radius of Earth, depth\n"
        inputd += "MFDECL  -27.05                            magnetic field declination (+E, -W)\n"

        spric = self.GetParameter('spric')
        if spric:
          par['spric'] = 'T'
        else:
          par['spric'] = 'F'

        if self.GetParameter('ranpri') == 3:
              inputd += "PNORM   %(pnormh)s %(pnormhe)s %(pnormn)s %(pnormal)s %(pnormfe)s      5-component relative contribution\n" 
              inputd += "PGAM    %(pgamh)s %(pgamhe)s %(pgamn)s %(pgamal)s %(pgamfe)s           5-component spectral indices\n" 

        self.steering += inputd
        
        # Write out DETPARAMS geometry file
        tmpdir    = par['tmpdir'] % par
        if not os.path.exists(tmpdir): os.makedirs(tmpdir)
        detparams = open(tmpdir+"/DETPARAMS",'w')
        detparams.write("-LENGTH=%(length)s -RADIUS=%(radius)s -DEPTH=%(depth)s\n" % par)
        detparams.close()


class ThinCorsika(Corsika):
    """
    This class provides an interface for preprocessing files in iceprod
    """

    def __init__(self):
        Corsika.__init__(self)
        self.name   = 'thincorsika'
        self.logger = logging.getLogger('iceprod::ThinCorsika')

        # ThinCorsika specific
        self.AddParameter('thinem_e','Fraction of primary energy where thinning algorithm is used for electromagnetic particles.',1.0E-6)
        self.AddParameter('thinem_wmax','maximum weight to be given to any thinned electromagnetic particle',10.0)  
        self.AddParameter('thinh_e','Energy(Gev) where thinning algorithm is used for hadrons',1.0)
        self.AddParameter('thinh_wmax','Maximum weight to be given to any thinned hadronic particle',1.0)
        
    def configure(self):
        Corsika.configure(self)

        par = self.parameters
        # Calculate parameters to be used for thinning
        par['efrcthn'] = par['thinem_e']
        par['wmax'] = par['thinem_wmax']
        par['rmax'] = 0.0
        par['thinrat'] = par['thinh_e']/par['thinem_e']
        par['weitrat'] = par['thinh_wmax']/par['thinem_wmax']
        
        inputd  = "THIN  %(efrcthn)f %(wmax)f %(rmax)f       EM thinning level weightmax rmax\n"
        inputd += "THINH  %(thinrat)f %(weitrat)f            Ratios for Hadronic thinning\n"

        self.steering += inputd


class AutoThinCorsika(Corsika):
    """
    This class provides an interface for preprocessing files in iceprod
    """

    def __init__(self):
        Corsika.__init__(self)
        self.name   = 'thincorsika'
        self.logger = logging.getLogger('iceprod::AutoThinCorsika')

        # ThinCorsika specific
        self.AddParameter('thin_method','Method for calculating thinning parameters.','2009')
        
    def configure(self):
        Corsika.configure(self)

        par = self.parameters
        if par['thin_method'] == '2009':
                # Calculate parameters to be used for thinning
                par['efrcthn'] = 1.0E-6
                par['wmax'] = par['emin']*par['efrcthn']    # Calculate optimum weight from Alessio
                if par['wmax'] < 1.0:   # Ensure max weight is at least 1
                   par['wmax'] = 1.0
                par['rmax'] = 0.0
                par['thinrat'] = 10.0/par['efrcthn']        # Just to be safe
                par['weitrat'] = 1.0/par['wmax']
        elif par['thin_method'] == '2010':
                # Calculate parameters to be used for thinning
                # No in-ice muon thinning

                par['efrcthn'] = 1.0E-6
                if par['emin'] > pow(10, 8.4) :
                   par['efrcthn'] = 273.0/par['emin']

                par['wmax'] = par['emin']*par['efrcthn']	

                if par['wmax'] < 1.0:	# Ensure max weight is at least 1
                    par['wmax'] = 1.0
                par['rmax'] = 0.0
                par['thinrat'] = 10.0/par['efrcthn']		# Just to be safe, ethem/ethhad
                par['weitrat'] = 1.0/par['wmax']                # wmaxem/wmaxhad
        else:
            self.logger.error('Specified thinning method not supported')

        inputd  = "THIN  %(efrcthn).3f %(wmax)f %(rmax)f       EM thinning level weightmax rmax\n"
        inputd += "THINH  %(thinrat)f %(weitrat)f            Ratios for Hadronic thinning\n"

        self.steering += inputd


class UCR(ipmodule.IPBaseClass):

   def __init__(self):
        ipmodule.IPBaseClass.__init__(self)
        self.logger = logging.getLogger('iceprod::UCR')

        self.AddParameter('ucr_binary','UCR executable','$I3_BUILD/bin/ucr-icetray-ucr')
        self.AddParameter('ucr_opts','UCR options','')
        self.AddParameter('input','UCR input','')
        self.AddParameter('output','UCR output','')


   def Execute(self,stats):
        if not ipmodule.IPBaseClass.Execute(self,stats): return 0

        ucr_bin   = self.GetParameter('ucr_binary')
        ucr_bin   = expandvars(ucr_bin)

        ucr_opts  = self.GetParameter('ucr_opts')
        ucr_opts  = expandvars(ucr_opts)

        ucr_in    = self.GetParameter('input')
        ucr_in    = expandvars(ucr_in)

        ucr_out   = self.GetParameter('output')
        ucr_out   = expandvars(ucr_out)

        ucr_cmd   = " ".join([ucr_bin,'-out='+ucr_out, ucr_in, ucr_opts])

        os.system("touch %s" % ucr_out) # for some reason ucr gets upset if the file doesn't exits
        # Run UCR 
        self.logger.info(ucr_cmd)
        status, output = getstatusoutput(ucr_cmd)

        self.logger.info(output)
        if status:
           self.logger.error(output)

        return status


