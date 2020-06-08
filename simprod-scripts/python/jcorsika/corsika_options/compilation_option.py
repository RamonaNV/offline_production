"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
 """

from option import OptionSuite, Option, OptionGroup
from . import coconut_options
from .get_version import get_version

class CompilationOption(OptionSuite):
    """
    FIXME : Write docs here.
    """
    def __init__(self, options, groups, version=None):
        OptionSuite.__init__(self, options)

        if version is None: version = get_version()

        from .coconut_options import CoconutOptionCollection
        self.coconut_options = CoconutOptionCollection(version)

        options.append(Option(long_name = '--version',
                              name = 'version',
                              help = 'Corsika version',
                              default = 'v6900'))

        options.append(Option(long_name = '--platform',
                              name = 'platform',
                              help = 'Compiler platform. This is usually just a tag for the tarball.',
                              default = 'Linux'))

        options.append(Option(long_name = '--fluka-tarball',
                              name = 'fluka_tarball',
                              help = 'Fluka tarball.',
                              default = None))

        options.append(Option(long_name = '--old-names',
                              name = 'old_names',
                              action = 'store_true',
                              help = 'Use old naming convention for executables\nUses corsika73700Linux_SIBYLL_gheisha instead of corsika73700Linux_1_6_1_1',
                              default = False))

        required = ('float', 'he_model', 'le_model', 'geometry')
        titles = {g:'' for g in self.coconut_options.option_groups}
        titles.update({ 'float':'32/64 bit compile option',
                        'he_model': 'High-energy hadronic model',
                        'le_model': 'Low-energy hadronic model',
                        'geometry': 'Detector geometry',
                        'other': 'Other compilation options'
                    })
        defaults = {g:None for g in self.coconut_options.option_groups}
        defaults.update({'float':'float_default'})
        for g in required:
            choices = [o.name for n,o in self.coconut_options.options.iteritems() if o.group==g]
            options.append(Option(g, long_name='--'+g.replace('_','-'),
                                  choices=choices,
                                  default=defaults[g],
                                  help=titles[g] + '. One of: %s'%(', '.join(choices))))
        for g in self.coconut_options.option_groups:
            if not g in required: groups.append(OptionGroup(name=g,
                                                            title=titles[g]))
        for n,o in self.coconut_options.options.iteritems():
            if not o.group in required:
                options.append(Option(long_name='--'+o.name.replace('_','-'),
                                      help=o.help,
                                      name=o.group,
                                      action='append_const',
                                      const=o.name,
                                      group=o.group))

    def binary_name(self):
        opt = self.coconut_options.get_coconut_options(self.options)
        if self.options['old_names']: name = self.coconut_options.old_name_from_options(opt)
        else: name = self.coconut_options.name_from_options(opt)
        return "corsika%sLinux_%s" % \
            (self.options['version'].lstrip('v'),
             #self.options['platform'],
             name
            )

    def configure(self, options):
        # call the base class method
        # this is not 'pure virtual' in the parlance of our times
        OptionSuite.configure(self, options)

        if self.options['other'] is None:
            self.options['other'] = []

    def validate(self):
        import logging
        res = True

        models = [o.name for o in self.coconut_options.options.values() if o.group=='he_model']
        valid = bool(self.options['he_model']) and self.options['he_model'].lower() in models
        if not valid: logging.warning('High-energy hadronic model not defined or not one of: %s'%(', '.join(models)))
        res *= valid

        models = [o.name for o in self.coconut_options.options.values() if o.group=='le_model']
        valid = bool(self.options['le_model']) and self.options['le_model'].lower() in models
        if not valid: logging.warning('Low-energy hadronic model not defined or not one of: %s'%(', '.join(models)))
        res *= valid
        if not res: logging.warn('invalid compilation option') 

        return res
