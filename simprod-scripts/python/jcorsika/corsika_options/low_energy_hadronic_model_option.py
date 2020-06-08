"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
 """

from option import OptionSuite, Option
import logging

class LowEnergyHadronicModelOption(OptionSuite):
    """
    FIXME : Write docs here.
    """
    def __init__(self, options, groups):
        OptionSuite.__init__(self, options)

        # this is set in the compilation options part
        #self.option_parser.add_option('--le_model',
        #                               dest = 'le_model',
        #                               help = 'Low-energy hadronic interaction model',
        #                               default = 'gheisha',
        #                               action = "store")

    def append_to_steering(self):
        return ''

    def validate(self):
        # main check is done in compilation_options
        import os
        # another hack... someone needs to set FLUPRO, but this requires corsika_path, which is set in configure, but because of deprecated_options it might be set after this suite, so it was bumped here...
        if self.options['le_model'] == 'fluka' and not 'FLUPRO' in os.environ:
            os.environ['FLUPRO'] = os.path.join(self.options['corsika_path'], 'fluka/')
        #if self.options['le_model'] == 'fluka' and not self.options['steering'] and not 'FLUPRO' in os.environ and not os.path.exists(self.options['corsika_path'] + 'fluka/'):
        #    logging.warn("Requested to use Fluka. You need to either define FLUPRO or make sure the directory fluka exists in corsika's location")
        #    return False
        return True
