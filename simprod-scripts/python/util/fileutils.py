import os
import subprocess
from .. import ipmodule

def download(url,output=None,username=None,password=None):
    """Download a file to output (default: basename of url)"""
    if not output:
        output = os.path.basename(url)
    if username and password:
        wget_cmd = 'wget -q --tries=2 --user=%s --password=%s --output-document=%%s %%s'%(username,password)
        curl_cmd = 'curl -s --retry 2 -u %s:%s -o %%s %%s'%(username,password)
    else:
        wget_cmd = 'wget -q --tries=2 --output-document=%s %s'
        curl_cmd = 'curl -s --retry 2 -o %s %s'
    globus_cmd = 'globus-url-copy -cd -r -nodcau -rst-retries 5 -rst-interval 60 %s %s'
    cp_cmd = 'cp %s %s'
    
    if url.startswith('http://') or url.startswith('ftp://'):
        if subprocess.call(wget_cmd%(output,url),shell=True):
            if subprocess.call(curl_cmd%(output,url),shell=True):
                raise Exception('cannot download %s'%url)
    elif url.startswith('gsiftp://'):
        if subprocess.call(globus_cmd%(url,output),shell=True):
            raise Exception('cannot download %s'%url)
    elif url.startswith('file:/'):
        url = url[5:]
        if output.startswith('file:/'):
            output = output[5:]
        if subprocess.call(cp_cmd%(url,output),shell=True):
            raise Exception('cannot download %s'%url)
    elif url.startswith('/'):
        if subprocess.call(cp_cmd%(url,output),shell=True):
            raise Exception('cannot download %s'%url)
    else:
        raise Exception('unknown download protocol for %s'%url)

def isurl(url):
    for prefix in ['http://','https://','ftp://','gsiftp://','file:/']:
        if url.startswith(prefix):
            return True
    return False

def untar(path):
    if path.endswith('.tgz') or path.endswith('.tar.gz'):
        return not subprocess.call('tar -zxf %s'%(path,),shell=True)
    elif path.endswith('.tar.bz2'):
        return not subprocess.call('tar -jxf %s'%(path,),shell=True)
    elif path.endswith('.tar'):
        return not subprocess.call('tar -xf %s'%(path,),shell=True)
    elif path.endswith('.zip'):
        return not subprocess.call('unzip %s'%(path,),shell=True)
    else:
        raise Exception('unknown archive format')


class SummaryMerger(ipmodule.ParsingModule):
    """
    This class provides an interface for preprocessing files in iceprod
    """

    def __init__(self):
        ipmodule.ParsingModule.__init__(self)
        self.AddParameter('inputfilelist','The names of the files you want to merge',[])
        self.AddParameter('outputfile','Name of output tarball',"summary.xml")


    def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0
        from . import ReadI3Summary, WriteI3Summary
        from icecube import dataclasses
        import math

        filelist = self.GetParameter('inputfilelist')
        outfile  = self.GetParameter('outputfile')
        retval = 0
        cmd = ''
        summary = dataclasses.I3MapStringDouble()
        for file in filelist:
            stats = ReadI3Summary(file)
            for key,value in stats.items():
                if isinstance(value,float) and math.isnan(value):
                   continue # don't add NaN
                if summary.has_key(key):
                    summary[key] += value
                else:
                    summary[key]  = value
        WriteI3Summary(summary, outfile)

        return 0

