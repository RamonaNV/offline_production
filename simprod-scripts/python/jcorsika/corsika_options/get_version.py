from optparse import OptionParser, BadOptionError
import logging

class _PassThroughOptionParser_(OptionParser):
    """
    A parser that just ignores options that throw exceptions when processed
    """
    def __init__(self, *args, **kwargs):
        self.quiet='quiet' in kwargs
        if 'quiet' in kwargs: del kwargs['quiet']
        OptionParser.__init__(self, *args, **kwargs)
    def _process_long_opt(self, rargs, values):
        try:
            OptionParser._process_long_opt(self, rargs, values)
        except Exception as e:
            if not self.quiet:
                logging.warn(str(e))
    def _process_short_opts(self, rargs, values):
        try:
            OptionParser._process_short_opts(self, rargs, values)
        except Exception as e:
            if not self.quiet:
                logging.warn(str(e))


class VersionError(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)

def get_version():
    from . import coconut_options
    available_keys = sorted(coconut_options.options_.keys())
    parser = _PassThroughOptionParser_('Arguments depend on corsika version. Do %s --version <version> to see specific options.', add_help_option=False, quiet=True)
    parser.add_option('--version', help='One of ' + ', '.join(available_keys))
    opts, args = parser.parse_args()
    if not opts.version:
        raise VersionError("CORSIKA version not specified. Use the --version option.\nUse one version of: %s\nExact options depend on the version. Use --version for more info."%(','.join(available_keys)))
    if not opts.version in available_keys:
        raise VersionError('CORSIKA version %s is not available. Only %s'%(opts.version, ', '.join(available_keys)))
    return opts.version
