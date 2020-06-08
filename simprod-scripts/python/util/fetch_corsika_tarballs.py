"""
 Utilities for fetching CORSIKA tarballs

 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: 
"""

# FIXME : Which imports aren't needed?

import os
import re
import sys
import math
import glob
import time
import string
import shutil
try:
    import cPickle as pickle
except:
    import pickle
from os import system
from os.path import expandvars
import getpass
import base64
import logging


# FIXME : What code is not meant for import and should be kept internal?
# Should only fetch_tarball be imported?

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
      logging.warning("sums differ")
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
        logging.info("retrieving %s..." % md5url)
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
            logging.info("retrieving %s..." % (meta['url'] + meta['suffix']))
            failed = wget(meta['url']+meta['suffix'])
            if failed:
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

