import os, sys
from os.path import expandvars
import time
import logging
import getpass
import base64

if sys.version_info[0] >= 3:
    import urllib.parse as urlparse
    from urllib.error import HTTPError
    from urllib.request import urlopen, Request
else:
    import urlparse
    from urllib2 import HTTPError
    from urllib2 import urlopen, Request

_logger = logging.getLogger('CorsikaStage')

def _strip_auth(url,auth=None):
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

def _wget(url, output_path, blocksize=2**16):
    auth = None
    _logger.info("fetching %s -> %s"%(url, output_path))

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
        f = urlopen(_strip_auth(url,auth=auth))
            
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

def _checksum(sum1,sum2):
   if not os.path.exists(sum1): return False
   if not os.path.exists(sum2): return False

   _logger.info("check md5sum %s %s"%(sum1, sum2))
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
      _logger.info("sums differ")
      return False
   else:
      _logger.info("sums match")
      return True

class _Lock:
    """
    Simple context-based lock
    """
    def __init__(self, filename):
        self.filename = filename
    def __enter__(self):
        try: # create lockfile
            os.mknod(self.filename)
        except OSError as oserr: # file exists
            _logger.warn("%s %s ." % (oserr, self.filename))
            time.sleep(300) # 5min
            try:
                if (time.time() - os.stat(self.filename).st_mtime) > 2400: 
                    os.remove(self.filename)
            except:
                pass
    def __exit__(self, *args):
        os.remove(self.filename)

class _TempCD:
    """
    Context manager to temporarily change to a directory.
    """
    def __init__(self, dirname):
        self.dirname = dirname
        self.back = os.getcwd()
    def __enter__(self):
        os.chdir(self.dirname)
    def __exit__(self, *args):
        os.chdir(self.back)

_tmpstorage = expandvars('$TMPDIR/etc/%s_icesoft' % getpass.getuser())
def fetch_tarball(url=None, filebase=None, suffix='.tar.gz', directory=_tmpstorage):
    """
    Get the corsika tarball

    The URL of the tarball is {url}/{filebase}{suffix}
    and its md5sum is {url}/{filebase}.md5sum
    This function assumes that the directory that results
    from uncompressing is called {filebase}

    Parameters
    ----------
    url: repository of tarballs
    filebase: the basename of the tarball without the suffix
    suffix: .tar.gz by default
    directory: the location where tarball is to be uncompressed

    Returns
    -------
    topdir: the directory that results from uncompressing.
    """
    meta = {'url':url, 'filebase':filebase, 'suffix':suffix}
    if ':' in directory:
        # more than one path specified, choose first good one
        storage = None
        for t in directory.split(':'):
            base = dirname(t)
            if not exists(base):
                continue
            if not exists(t):
                try: os.makedirs(t)
                except: continue
            storage = t
            break
        if storage:
            directory = storage
        else:
            directory = _tmpstorage

    lockfile = os.path.join(directory,"%s.lock" % meta['filebase'])

    cachedir = os.path.join(directory,meta['filebase'])

    md5sum_src = meta['url'] + '.md5sum'
    md5sum_ref = os.path.join(directory, meta['filebase'] + '.md5sum.ref')
    md5sum_new = os.path.join(directory, meta['filebase'] + '.md5sum')
    tarball_src = meta['url'] + meta['suffix']
    tarball_dest = os.path.join(directory, meta['filebase'] + meta['suffix'])

    _logger.debug("file sources:\n  %s\n  %s"%(tarball_src, md5sum_src))
    _logger.debug("file destinations:\n  %s\n  %s"%(tarball_dest, md5sum_ref))

    #fetch checksum
    _wget(md5sum_src, md5sum_ref)
    # See if we need to download new tarball
    trials = 5 # limit number of attempts
    while  (not os.path.exists(tarball_dest)) or \
           (not os.path.exists(cachedir)) or \
           (not _checksum(md5sum_ref, md5sum_new)): 
        _logger.debug("trials left: %d"%trials)
        with _Lock(lockfile):
            if os.system("rm -rf %s %s" % (tarball_dest,cachedir)):
                os.system('/usr/bin/printenv')
                raise Exception('failed to remove cached %s'%tarball_dest)

            failed = _wget(tarball_src, tarball_dest)

            if not os.path.isfile(failed):
                raise Exception("unable to fetch file from %s" % tarball_src)

            if os.path.isdir(meta['filebase']): # clean old i3_work
                if os.system("rm -rf " + meta['filebase']):
                    raise Exception('failed to remove %s'%(meta['filebase']))

            # create md5sum
            if os.system("tar -C%s -xzf %s" % (directory, tarball_dest)):
                raise Exception("failed to untar %s"%tarball_dest)
            with _TempCD(os.path.dirname(tarball_dest)):
                # this chdir is necessary because of the simple way how md5sum is checked.
                if os.system("md5sum %s > %s" % (os.path.basename(tarball_dest), md5sum_new)):
                    raise Exception("failed to calculate md5sum of %s"%(tarball_dest))

            trials -= 1 # decrement counter
            if trials < 1: 
                raise Exception("Max no. of trails reached to stage metaproject '%s'" % meta['filebase'])

    return os.path.join(directory,meta['filebase'])

def corsika_stage(version, url=None, location=None, platform=None, directory=os.environ['PWD'], suffix='.tar.gz'):
    """
    Prepare CORSIKA files and executables.

    As a result of this function, there should be a corsika directory
    containing the run or bin directory where CORSIKA executables can
    be found. This corsika directory is expected to be named according
    to this formats:
        corsika-{version}
        corsika-{version}.{platform}
    This function first looks for the corsika directory in the location
    specified by the 'location' parameter. If not found, it then looks
    for a tarball at the following URLs:
        {url}/corsika-{version}{suffix}
        {url}/corsika-{version}.{platform}{suffix}

    Parameters
    ----------
    version: corsika version (something like '73700')
    url: A url where corsika tarballs can be downloaded.
    location: A directory where one should find the corsika directory.
    platform: This is normally just an identifying string added to the file or directory name.
    directory: The directory where the corsika directory should end up.
    suffix: The suffix of the tarball. '.tar.gz' by default. This is replaced by .md5sum to find the corresponding md5sum file.

    Returns
    -------
    topdir : string. The path to the corsika location (the directory containing the run or bin directory)
    """

    directory = os.path.abspath(directory)
    if not os.path.exists(directory):
        os.makedirs(directory)

    if platform: platform = '.'+platform
    filebase = "corsika-%s%s" % (version, platform)
    url = os.path.join(url,filebase)

    if location:
        location_path = os.path.join(location,filebase)
        _logger.info("Corsika location configured to '%s'",location)

    if location and os.path.exists(location_path):
        _logger.info("Corsika location found: '%s'",location_path)
        topdir = location_path
    else:
        if location: _logger.warn('Corsika location not found (%s), will try to download.', location_path)
        _logger.info("Downloading tarball: '%s'", url + suffix)
        topdir = fetch_tarball(url=url, filebase=filebase, suffix=suffix, directory=directory)

    return topdir

        
if __name__ == "__main__":
    #test:
    logging.basicConfig(level=logging.INFO)
    import corsika_dev.corsika_stage
    r = corsika_dev.corsika_stage.corsika_stage(url='file:/data/sim/sim-new/downloads',
                                                version='v73700',
                                                platform='gfortran_4.8.2_IT')
    print(r)

