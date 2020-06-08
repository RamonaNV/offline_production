"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
"""

import logging
import subprocess
import os, sys
from os.path import expandvars

class CorsikaBinary:
    """
    In the simplest case a CorsikaBinary should only have take a
    binary name and execute that with the given INPUTS file.  It
    might need to do other stuff to work around CORSIKA's pathology,
    but we'll see...
    This only requires the bare minimun...you have a INPUTS file
    and a binary, so go!
    We might want to add back in checks, cleanups, etc...
    """
    def __init__(self, options, binary_name):
        """
        I don't know what this is supposed to do yet.
        """
        self.options = options
        self.logger = logging.getLogger('Corsika')
        # staging goes here
        # When we make tarballs we put the executables in a 'bin' directory
        # The default for straight unadulterated CORSIKA builds are in 'run'
        # So we check both.
        if os.path.exists(os.path.join(self.options['corsika_path'], 'bin/')) :
            self.executable_path = os.path.join(self.options['corsika_path'], 'bin/')
        elif os.path.exists(os.path.join(self.options['corsika_path'], 'run/')) :
            self.executable_path = os.path.join(self.options['corsika_path'], 'run/')
        else:
            logging.critical("Executables are expected to be in either a 'run' or 'bin' path.")
            logging.critical("Can't find either in '%s'" % self.options['corsika_path'])
            raise Exception

        self.absolute_path = os.path.abspath(self.executable_path + binary_name)



    def Execute(self, steering_filename):
        # Before we execute, let's make sure all the pieces are in-place
        # and we're setting Popen up to succeed.

        if os.path.exists(self.options['corsika_path'] + 'lib/') :
            os.environ['LD_LIBRARY_PATH'] = self.options['corsika_path'] + 'lib:' + os.environ['LD_LIBRARY_PATH']
        if not os.path.exists(steering_filename):
            logging.critical("Steering file '%s' does not exist." % steering_filename)
            raise Exception
        if not os.path.exists(self.options['outdir']):
            logging.critical("outdir dir '%s' does not exist." % self.options['outdir'])
            raise Exception
        if not os.path.exists(self.absolute_path):
            logging.critical("Binary file '%s' does not exist." % self.absolute_path)
            raise Exception

        # Run CORSIKA
        steering_file = open(steering_filename)
        log_file = None
        if self.options['log']: log_file = open(self.options['log_filename'],'w')
        self.logger.info('%s < %s', self.absolute_path, steering_filename)
        corsika_process = subprocess.Popen(self.absolute_path,
                                           stdin = steering_file,
                                           stdout = log_file if log_file else sys.stdout,
                                           stderr = log_file if log_file else sys.stderr,
                                           cwd = self.executable_path)
        if log_file: log_file.close()

        output = corsika_process.communicate()[0]
        steering_file.close()

        return 0

