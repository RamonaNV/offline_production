#!/usr/bin/env python
import logging
logging.basicConfig()
rootLogger = logging.getLogger('')
rootLogger.setLevel(logging.INFO)

from icecube.simprod.modules import Corsika5ComponentGenerator

if __name__ == '__main__':
   stats = {}
   cors = Corsika5ComponentGenerator()
   cors.ExecuteOpts(stats)
