from .low_energy_hadronic_model_option import LowEnergyHadronicModelOption
from .high_energy_hadronic_model_option import HighEnergyHadronicModelOption
from .execution_option import ExecutionOption
from .compilation_option import CompilationOption
from .parameters_option import ParametersOption
from .atmosphere_option import AtmosphereOption
from .thinning_option import ThinningOption
from .deprecated_options import DeprecatedOptions
import get_version
import logging

class CorsikaOptions:
    """
    This class manages an OptionParser and a collection
    of Options, which determine the configuration. Each
    Option is responsible for it's own sector of the
    configuration of CORSIKA.

    This class creates and binary and executes it creating
    an output file.
    """
    def __init__(self, option_parser=None, defaults={}, version=None):
        self.logger = logging.getLogger('CorsikaOptions')
        self._option_parser = option_parser
        self.options = {}
        options = []
        groups = []
        try:
            self.corsika_options = [
                ExecutionOption(options, groups),
                CompilationOption(options, groups, version),
                ParametersOption(options, groups),
                HighEnergyHadronicModelOption(options, groups),
                LowEnergyHadronicModelOption(options, groups),
                AtmosphereOption(options, groups),
                ThinningOption(options, groups),
                DeprecatedOptions(options, groups)
            ]
        except get_version.VersionError as e:
            self.logger.critical("Exception: "+str(e))
            raise e
        for i in range(len(options)):
            if options[i].name in defaults:
                options[i].default = defaults[options[i].name]
        self._optparse_options = options
        self._optparse_groups = groups
        if self._option_parser: self.parse_command_line()
        else: self.update({o.name:o.default for o in options})

    def binary_name(self):
        return self.corsika_options[1].binary_name()

    def parse_command_line(self, parser=None):
        from optparse import OptionGroup
        from .get_version import _PassThroughOptionParser_ as OptionParser
        if parser is None: self._option_parser = OptionParser()
        option_groups = {g.name:OptionGroup(self._option_parser, g.title, g.help) for g in self._optparse_groups}
        for o in self._optparse_options:
            args = [n for n in [o.short_name, o.long_name] if n]
            kwargs = {'dest':o.name,
                      'help':o.help,
                      'default':o.default,
                      'action':o.action,
                      'const':o.const,
                      'choices':o.choices,
                      'type': o.type}
            kwargs = dict([(k,v) for k,v in kwargs.iteritems() if not v is None])
            if o.group:
                option_groups[o.group].add_option(*args, **kwargs)
            else:                
                self._option_parser.add_option(*args, **kwargs)
        for n,g in option_groups.iteritems():
            self._option_parser.add_option_group(g)

        (options, args) = self._option_parser.parse_args()
        self.update(vars(options))

    def update(self, options):
        if not self.options: self.options = {}
        self.options.update(options)
        self.logger.debug(str(options))
        for option in self.corsika_options:
            option.configure(self.options)

    def is_valid(self):
        return not bool(sum([(not option.validate()) for option in self.corsika_options]))

    def steering(self):
        steering_str = ""
        for option in self.corsika_options:
            steering_str += option.append_to_steering()
        steering_str += "EXIT \n"
        return steering_str

    def __len__(self):
        if not self.options: return 0
        return len(self.options)

    def __getitem__(self, k):
        return self.options[k]

    def __contains__(self, k):
        return bool(self.options) and k in self.options

    def __iter__(self):
        return self.options.__iter__()
