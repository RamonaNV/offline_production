"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
"""
import os
import sys
import logging

from .corsika_binary import CorsikaBinary
from .corsika_builder import CorsikaBuilder
from .corsika_options import CorsikaOptions

class Corsika:
    """
    This class contains a CorsikaOptions, creates and binary and executes it creating
    an output file.
    """
    def __init__(self, option_parser=None, defaults={}, version=None):
        self.config = CorsikaOptions(option_parser, defaults=defaults, version=version)

    def WriteSteering(self):
        """
        Writes a steering file
        """
        if not os.path.exists(self.config['outdir']):
            os.makedirs(self.config['outdir'])
        self.steering_filename = self.config['outdir'] + "DAT%06d.in" % self.config['runnr']
        inputfile = open(self.steering_filename,'w')
        inputfile.write(self.config.steering())
        inputfile.close()

    def Execute(self, stats={}):
        """
        Run CORSIKA
        """
        if not self.config.is_valid():
            raise Exception("Invalid configuration")

        if self.config['steering']:
            print(self.config.steering())
            return

        binary_name = self.config.binary_name()
        if len(binary_name) == 0 :
            logging.critical("Empty binary name.")
            raise Exception("Empty binary name.")

        if not os.path.exists(self.config['corsika_path']):
            logging.info('CORSIKA path does not exist. About to fetch corsika tarball')
            import corsika_dev.corsika_stage
            location = corsika_dev.corsika_stage.corsika_stage(url=self.config['url'],
                                                               version=self.config['version'],
                                                               platform=self.config['platform'])
            if location != self.config['corsika_path']:
                logging.warn("Inconsistent CORSIKA paths. The configured value ({0}) is not the same as the downloaded one ({1}). Weird things might happen.".format(self.config['corsika_path'], location))

        # unfortunately, a hack. One cannot know whether the directory is there until after staging:
        if self.config['le_model'] == 'fluka' and (not 'FLUPRO' in os.environ or not os.path.exists(os.environ['FLUPRO'])):
            logging.warn("Requested to use Fluka. You need to either define FLUPRO or make sure the directory fluka exists in corsika's location")
            raise Exception("FLUPRO undefined or non-existent (%s)"%(os.environ['FLUPRO'] if 'FLUPRO' in os.environ else None))

        # Now simply instantiate a CorsikaBinary object and run
        binary = CorsikaBinary(self.config, binary_name)

        if not os.path.exists(binary.absolute_path) and self.config['build']:
            self.Build()

        if not os.path.exists(binary.absolute_path):
            logging.critical("Binary file '%s' does not exist." % binary.absolute_path)
            import glob
            v = self.config['version']
            for f in glob.glob(os.path.join(binary.executable_path, 'corsika%s*'%v)):
                logging.critical("  Available: %s"%(f))
            raise Exception

        self.WriteSteering()
        binary.Execute(self.steering_filename)


    def Build(self):
        """
        Build CORSIKA
        """
        builder = CorsikaBuilder(self.config['version'])
        filename = os.path.join(self.config['corsika_path'], "corsika_build.log")
        if not self.config['dryrun']:
            log_file = open(filename, "r+" if os.path.exists(filename) else "w")
            log_file.write("corsika.py building...")

        if self.config['le_model']=='fluka':
            fluka_dir = os.path.join(self.config['corsika_path'], 'fluka')
            if not os.path.exists(fluka_dir):
                logging.info("Building Fluka")
                if not self.config['fluka_tarball']:
                    raise Exception('If you want Fluka, you have to specify the --fluka-tarball option.')
                if not self.config['dryrun']: builder.build_fluka(fluka_dir, self.config['fluka_tarball'], log_file=log_file)
            os.environ['FLUPRO'] = fluka_dir
            os.environ['F77']='gfortran'
        logging.info("Building CORSIKA")
        if not self.config['dryrun']:
            builder.coconut(builder.options.get_coconut_options(vars(self.config)), log_file, cwd=self.config['corsika_path'])
            log_file.close()

