"""
 copyright  (c) 2017 the icecube collaboration

 A collection of IceProd modules to run CORSIKA.
 There is one module for each supported version.

 @version: $Revision: $
 @date: $Date: $
 @author: Javier Gonzalez <javierg@udel.edu>
"""
from icecube.simprod import ipmodule

class _CorsikaSimprod(ipmodule.IPBaseClass):
    def __init__(self, version):
        from .corsika import Corsika
        ipmodule.IPBaseClass.__init__(self)
        self.corsika = Corsika(version=version)
        multiple = ['other']
        for o in self.corsika.config._optparse_options:
            if o.name in multiple:
                ipmodule.IPBaseClass.AddParameter(self, o.const, '', False)
            else:
                ipmodule.IPBaseClass.AddParameter(self, o.name, '', o.const)
    def Execute(self, options, stats={}):
        self.logger.info("executing")
        if not ipmodule.IPBaseClass.Execute(self,stats): return 0
        multiple = ['other']
        options = {k:[] for k in multiple}
        for o in self.corsika.config._optparse_options:
            if o.name in multiple:
                self.logger.debug(' %20s: %s', o.const, bool(self.GetParameter(o.const)))
                if self.GetParameter(o.const): options[o.name].append(o.const)
            else:
                self.logger.debug(' %20s %s', o.name, self.GetParameter(o.name))
                if not self.GetParameter(o.name) is None: options[o.name] = self.GetParameter(o.name)
        self.corsika.config.update(options)
        return self.corsika.Execute(stats)
    def Finish(self):
        self.logger.info("finishing")

_class_def = """
class Corsika{version}(_CorsikaSimprod):
    def __init__(self):
        _CorsikaSimprod.__init__(self, version='{version}')"""

from .corsika_options import coconut_options
for _k in coconut_options.options_.keys():
    exec(_class_def.format(version=_k))
del coconut_options
