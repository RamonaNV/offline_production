"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Javier Gonzalez <javierg@udel.edu>
 """

from option import OptionSuite, Option, OptionGroup

class DeprecatedOptions(OptionSuite):
    """
    A container for deprecated options.
    """
    def __init__(self, options, groups):
        OptionSuite.__init__(self, options)

        groups.append(OptionGroup(name='deprecated',
                                  title="Deprecated options",
                                  help="Stop using these!"))

        options.append(Option(long_name='--runnum',
                              help = 'Run number. Use runnr instead.', group='deprecated'))
        options.append(Option(long_name='--crtype',
                              help = 'Primary particle code. Use primpar instead.', group='deprecated'))
        options.append(Option(long_name='--model', help = 'Use he_model instead.', group='deprecated'))
        options.append(Option(long_name='--lemodel', help = 'Use le_model instead.', group='deprecated'))
        options.append(Option(long_name='--nevents', help = 'Use nshow instead.', group='deprecated'))
        options.append(Option(long_name='--flat-detector', help = 'Use geometry instead.', group='deprecated'))
        options.append(Option(long_name='--usepipe', help = 'Use pipe instead.', group='deprecated'))
        options.append(Option(long_name='--topdir', help = 'Use outdir instead.', group='deprecated'))
        options.append(Option(long_name='--URL', help = 'Use url instead.', group='deprecated'))


    def configure(self, options):
        # call the base class method
        # this is not 'pure virtual' in the parlance of our times
        OptionSuite.configure(self, options)

        deprecated = []
        if not self.options['runnum'] is None:
            deprecated.append('runnum')
            self.options['runnr'] = self.options['runnum']
        if not self.options['crtype'] is None:
            deprecated.append('crtype')
            self.options['prmpar'] = self.options['crtype']
        if not self.options['model'] is None:
            deprecated.append('model')
            self.options['he_model'] = self.options['model']
        if not self.options['lemodel'] is None:
            deprecated.append('lemodel')
            self.options['le_model'] = self.options['lemodel']
        if not self.options['nevents'] is None:
            deprecated.append('nevents')
            self.options['nshow'] = self.options['nevents']
        if not self.options['flat_detector'] is None:
            deprecated.append('flat_detector')
            self.options['geometry'] = 'flat'
        if not self.options['usepipe'] is None:
            deprecated.append('usepipe')
            self.options['pipe'] = self.options['usepipe']
        if not self.options['topdir'] is None:
            deprecated.append('topdir')
            self.options['outdir'] = self.options['topdir']
        if not self.options['URL'] is None:
            deprecated.append('URL')
            self.options['url'] = self.options['URL']
        if deprecated:
            import logging
            logging.getLogger('Corsika').warning('Update your configuration. The following deprecated options are being used:\n  %s'%('\n  '.join(deprecated)))
