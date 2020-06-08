
import sys
import os
import math
from .. import ipmodule


class Splitter(ipmodule.IPBaseClass):
    """Split an input file into multiple parts.
       Round-robin splitting until all frames are output."""
    def __init__(self):

        ipmodule.IPBaseClass.__init__(self)
        self.AddParameter('inputfile','input i3 file','')
        self.AddParameter('outfileprefix','output file prefix','')
        self.AddParameter('noutputs','number of output files to split into',10)
        self.AddParameter('outfilesuffix','output file suffix (file extension)','.i3.bz2')

    def Execute(self,stats):
        if not ipmodule.IPBaseClass.Execute(self,stats): return 0
        
        from icecube import icetray, dataclasses, dataio
        
        infile  = self.GetParameter('inputfile')
        if not os.path.isfile(infile):
            raise Exception('input file does not exist')
        prefix = self.GetParameter('outfileprefix')
        noutputs = self.GetParameter('noutputs')
        if noutputs <= 1:
            raise Exception('noutputs must be greater than 1')
        suffix = self.GetParameter('outfilesuffix')
        
        # make output files
        outfiles = []
        pattern = '%s%0'+('%d'%math.ceil(math.log10(noutputs)))+'d%s'
        for i in range(noutputs):
            outfiles.append(dataio.I3File(pattern%(prefix,i,suffix),'w'))

        # read from file
        for n,frame in enumerate(dataio.I3File(infile)):
            outfiles[n%noutputs].push(frame)
        
        # close output files
        for file in outfiles:
            file.close()
        
        return 0


class Combiner(ipmodule.IPBaseClass):
    """Combine multiple input files into one output file.
       Grabs 1 frame from each file round-robin until all files
       are exhausted.  Designed as a compliment to the Splitter."""
    def __init__(self):

        ipmodule.IPBaseClass.__init__(self)
        self.AddParameter('outputfile','output i3 file','')
        self.AddParameter('infileprefix','input file prefix','')
        self.AddParameter('ninputs','number of input files to merge',10)
        self.AddParameter('infilesuffix','input file suffix (file extension)','.i3.bz2')

    def Execute(self,stats):
        if not ipmodule.IPBaseClass.Execute(self,stats): return 0
        
        from icecube import icetray, dataclasses, dataio
        
        outfile  = self.GetParameter('outputfile')
        prefix = self.GetParameter('infileprefix')
        ninputs = self.GetParameter('ninputs')
        if ninputs <= 1:
            raise Exception('ninputs must be greater than 1')
        suffix = self.GetParameter('infilesuffix')
        
        # make list of input files
        infiles = []
        pattern = '%s%0'+('%d'%math.ceil(math.log10(ninputs)))+'d%s'
        for i in range(ninputs):
            finput = pattern%(prefix,i,suffix)
            if not os.path.isfile(finput):
                raise Exception('file %s not found'%finput)
            infiles.append(dataio.I3File(finput))

        # read from files
        file = dataio.I3File(outfile,'w')
        nn = range(0,ninputs)
        try:
            while infiles:
                # try to read another frame off each file
                deln = []
                for n,infile in enumerate(infiles):
                    frame = None
                    try:
                        frame = infile.pop_frame()
                    except:
                        # must be out of frames
                        infile.close()
                        deln.append(n)
                    else:
                        if not frame:
                            # must be out of frames
                            infile.close()
                            deln.append(n)
                        else:
                            file.push(frame)
                for n in deln:
                    del infiles[n]
        except Exception as e:
            icetray.logging.log_error('Exception when combining')
            raise
        
        # close output file
        file.close()
        
        return 0
