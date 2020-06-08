#!/usr/bin/env python
import logging
logging.basicConfig()
rootLogger = logging.getLogger('')
rootLogger.setLevel(logging.INFO)

from icecube.simprod.modules import PolyplopiaMCPEMerge

if __name__ == '__main__':
   stats = {}
   merge = PolyplopiaMCPEMerge()
   merge.ExecuteOpts(stats)
