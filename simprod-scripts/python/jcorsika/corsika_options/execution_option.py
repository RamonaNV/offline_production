"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
"""
import os
import getpass
import logging

from option import OptionSuite, Option, OptionGroup

class ExecutionOption(OptionSuite):
    """
    Execution options that have nothing to do with CORSIKA itself (paths, logging, testing)
    """
    def __init__(self, options, groups):
        OptionSuite.__init__(self, options)

        groups.append(OptionGroup(name='execution',
                                  title='Execution Options',
                                  help='These are options dealing with the process itself'))

        options.append(Option(long_name='--outfile',
                              name = 'outfile',
                              help = 'Output file',
                              default = None,
                              group='execution'))

        options.append(Option(long_name='--log',
                              name = 'log',
                              help = 'Pipe stdout and stderr to file',
                              default = False,
                              action = "store_true",
                              group='execution'))

        options.append(Option(long_name='--compress',
                              name = 'compress',
                              help = 'Compress the output file.',
                              default = False,
                              action = "store_true",
                              group='execution'))

        options.append(Option(long_name='--build',
                              name = 'build',
                              help = 'Build executable if not present',
                              default = False,
                              action = "store_true",
                              group='execution'))

        options.append(Option(long_name='--pipe',
                              name = 'pipe',
                              help = 'write to named pipe',
                              default = False,
                              action = "store_true",
                              group='execution'))

        options.append(Option(long_name='--outdir',
                              name = 'outdir',
                              help = 'output directory',
                              default = os.environ['PWD'],
                              group='execution'))

        options.append(Option(long_name='--corsika-path',
                              name = 'corsika_path',
                              help = 'Path to corsika binary, files, etc... (this directory must contain the run or bin directory)', 
                              default = None,
                              group='execution'))

        options.append(Option(long_name='--url',
                              name = 'url',
                              help = 'fetch tarball from URL',
                              default = '',
                              group='execution'))

        options.append(Option(long_name='--cache',
                              name = 'cache',
                              help = 'Should cache tarball?',
                              default = False,
                              action = "store_true",
                              group='execution'))

        options.append(Option(long_name='--dryrun',
                              name = 'dryrun',
                              help = 'Only generate INPUTS file and exit',
                              default = False,
                              action = "store_true",
                              group='execution'))

        options.append(Option(long_name='--steering',
                              name = 'steering',
                              help = 'Write input on the screen',
                              default = False,
                              action = "store_true",
                              group='execution'))

        options.append(Option(long_name='--verbose',
                              name = 'verbosity',
                              help = 'Write INFO logging messages',
                              default = logging.WARNING,
                              action = "store_const",
                              const = logging.INFO,
                              group='execution'))


    def configure(self, options):
        # call the base class method
        # this is not 'pure virtual' in the parlance of our times
        OptionSuite.configure(self, options)

        logging.getLogger('Corsika').setLevel(options['verbosity'])

        # Make sure all directories end with a trailing '/'
        for attribute in ('outdir', 'corsika_path'):
            if self.options[attribute]:
                self.options[attribute] = os.path.abspath(self.options[attribute])
            if self.options[attribute] and \
              not self.options[attribute].endswith('/'):
                self.options[attribute] = self.options[attribute] + '/'

        # the default value of the outfile name depends on the runnum,
        # which is not a member of this policy, so it has to be set here.
        if not self.options['outfile']:
            self.options['outfile'] = "DAT%06d" % self.options['runnr']
        if self.options['log']:
            self.options['log_filename'] = self.options['outdir'] + "/DAT%06d.log" % self.options['runnr']
            logging.getLogger('Corsika').warn('log file: %s'%self.options['log_filename'])

    def validate(self):
        valid = True
        # another hack. Things should not be set here, but URL is a deprecated option that could be defined later
        if not self.options['corsika_path']:
            if self.options['url']:
                # will download
                self.options['corsika_path'] = '{pwd}/corsika-{version}.{platform}/'.format(pwd=os.environ['PWD'], **self.options)
            elif not self.options['steering']:
                msg = "You need to either specify the path to CORSIKA or a URL to fetch a tarball."
                valid = False
        return valid

    def append_to_steering(self):

        steering_str = ""

        if self.options['pipe'] :
            steering_str += "PIPE T \n"

        steering_str += "DIRECT %s \n" % self.options['outdir']

        return steering_str
