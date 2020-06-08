#!/usr/bin/env python
import logging
logging.basicConfig()
rootLogger = logging.getLogger('')
rootLogger.setLevel(logging.INFO)

from icecube.simprod.modules import CorsikaGenerator
from icecube.simprod.modules import Corsika5ComponentGenerator

if __name__ == '__main__':
   stats = {}
   cors = CorsikaGenerator()
   cors.ExecuteOpts(stats)
   print(stats)
