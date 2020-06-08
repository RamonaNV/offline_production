"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Alex Olivas <olivas@icecube.umd.edu>
"""

from .. import ipmodule
import logging
import subprocess

from .corsika_binary import CorsikaBinary

class CorsikaBinaryStager(CorsikaBinary):
    """
    This module manages a configuration and a corsika binary
    """
    def __init__(self, configuration):
        CorsikaBinary.__init__(self, configuration)

    def _Stage(self):        
        # FIXME : move all the staging code here
        # need to set self.binary_name and self.binary_path
        # accordingly.

        # First however, move and document the staging code
        # to a location in...say utils with a name like
        # corsika_staging_utils        
        pass
    
    def Execute(self):
        self._Stage() # this sets self.binary_name and self.binary_path   
        CorsikaBinary.Execute(self)
