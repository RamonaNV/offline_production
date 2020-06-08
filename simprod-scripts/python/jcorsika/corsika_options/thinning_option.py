"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
 """

from option import OptionSuite, Option, OptionGroup
import logging

class ThinningOption(OptionSuite):
    """
    Option to set thinning parameters in CORSIKA.

    This module replicates the behavior of the older simprod ThinCorsika and AutoThinCorsika.
    """
    def __init__(self, options, groups):
        OptionSuite.__init__(self, options)

        # this correspond one-to-one with corsika keywords
        options.append(Option(long_name = '--thin-fraction',
                              default = 1e-6,
                              help = 'Relative thinning factor (E_thin/E_0)'))
        options.append(Option(long_name = '--thin-wmax',
                              default = None,
                              help = 'Maximum weight'))
        options.append(Option(long_name = '--thin-rmax',
                              default = 0.,
                              help = 'Maximum weight. Default is 0 so the core is not thinned.'))
        options.append(Option(long_name = '--thin-thinrat',
                              default = 1.,
                              help = 'Ratio between thinning factor of EM and hadron components set using THINH. Default is 1.'))
        options.append(Option(long_name = '--thin-weitrat',
                              default = None,
                              help = 'Ratio between maximum weights of EM and hadron components set using THINH. The default depends on the method.'))

        # these are options that were in ThinCorsika
        # one problem with this is that the parameters do not correspond to CORSIKA keywords,
        # another is that the defaults are absurd.
        # Normally one does not need to give both (EM, hadron) thinning factors.
        options.append(Option(long_name = '--thinem_e',
                              help = 'Relative thinning energy for EM particles (deprecated)'))
        options.append(Option(long_name = '--thinem_wmax',
                              help = 'Maximum weight for EM particles (deprecated)'))
        options.append(Option(long_name = '--thinh_e',
                              help = 'Relative thinning energy for hadronic component (default is same as for EM) (deprecated)'))
        options.append(Option(long_name = '--thinh_wmax',
                              help='Maximum weight for particles in the hadronic component (default is 1) (deprecated)'))

        # these were options in AutoThinCorsika
        options.append(Option(long_name = '--thin_method',
                              help = 'Method for calculating thinning parameters. One of: none, 2009, 2010',
                              choices = [ 'none', '2009', '2010', 'auger-optimum', 'icecube-optimum'],
                              default = 'auger-optimum'))


    def configure(self, options):
        # call the base class method
        # this is not 'pure virtual' in the parlance of our times
        OptionSuite.configure(self, options)

        if not 'thin' in self.options['other']:
            return

        deprecated = [k for k in ['thinem_e', 'thinem_wmax', 'thinh_e', 'thinh_wmax'] if not self.options[k] is None]
        if deprecated:
            logging.getLogger('Corsika').warning('Update your configuration. The following deprecated thinning options are being used:\n  %s'%('\n  '.join(deprecated)))

        if not self.options['thinem_e'] is None:
            if self.options['thin_fraction'] is None:
                self.options['thin_fraction'] = self.options['thinem_e']
            else:
                logging.getLogger('Corsika').warning('Ignoring deprecated option thinem_e. Using thin-fraction instead.')

        if not self.options['thinh_e'] is None or not self.options['thinem_e'] is None:
            self.options['thin_thinrat'] = self.options['thinh_e']/self.options['thinem_e']
        if not self.options['thinh_wmax'] is None or not self.options['thinem_wmax'] is None:
            self.options['thin_weitrat'] = self.options['thinh_wmax']/self.options['thinem_wmax']

        if str(self.options['thin_method']) == 'none':
            self.options['efrcthn'] = self.options['thin_fraction']
            self.options['wmax'] = self.options['thin_wmax']
            self.options['rmax'] = self.options['thin_rmax']
            self.options['thinrat'] = self.options['thin_thinrat']
            self.options['weitrat'] = self.options['thin_weitrat']
        elif self.options['thin_method'] == '2009':
            self.options['efrcthn'] = 1.0E-6
            self.options['wmax'] = self.options['emin']*self.options['efrcthn'] # Calculate optimum weight from Alessio
            if self.options['wmax'] < 1.0: self.options['wmax'] = 1.0        # Ensure max weight is at least 1
            self.options['rmax'] = 0.0
            self.options['thinrat'] = 10.0/self.options['efrcthn']           # Just to be safe
            self.options['weitrat'] = 1.0/self.options['wmax']
        elif self.options['thin_method'] == '2010':
            # No in-ice muon thinning
            self.options['efrcthn'] = 1.0E-6
            # this is ~ E/(efrcthn * 10^8.4). Thinning is then 'absolute' with Eth=273 GeV if E>10**8.4 GeV
            if self.options['emin'] > pow(10, 8.4): self.options['efrcthn'] = 273.0/self.options['emin']
            self.options['wmax'] = self.options['emin']*self.options['efrcthn']
            if self.options['wmax'] < 1.0: self.options['wmax'] = 1.0  # Ensure max weight is at least 1
            self.options['rmax'] = 0.0
            self.options['thinrat'] = 10.0/self.options['efrcthn']     # Just to be safe, ethem/ethhad
            self.options['weitrat'] = 1.0/self.options['wmax']         # wmaxem/wmaxhad
        elif self.options['thin_method'] == 'auger-optimum':
            # M. Kobal (P. Auger Collaboration), Astropart. Phys. 15 (2001) 259
            # WMAX = EFRCTHN * Eo/GeV, WMAX >= 1
            self.options['efrcthn'] = self.options['thin_fraction']
            self.options['wmax'] = self.options['emin']*self.options['efrcthn']
            if self.options['wmax'] < 1.0: self.options['wmax'] = 1.0  # Ensure max weight is at least 1
            self.options['rmax'] = self.options['thin_rmax']
            self.options['thinrat'] = self.options['thin_thinrat']
            if self.options['thin_weitrat'] is None:
                self.options['weitrat'] = 1./100.
            else:
                self.options['weitrat'] = self.options['thin_weitrat']
        elif self.options['thin_method'] == 'icecube-optimum':
            # M. Kobal (P. Auger Collaboration), Astropart. Phys. 15 (2001) 259
            # WMAX = EFRCTHN * Eo/GeV, WMAX >= 1 and hadronic part is not thinned
            # WEITRAT = WMAX (with this, one can still get weights like 1.01 or so)
            # THINRAT > EFRCTHN * Eo / Ecut, where Ecut is low enough so muons/pions are never thinned. We take the ectual cut used.
            ecut = min(self.options['ecuts1'], self.options['ecuts2'])
            self.options['efrcthn'] = self.options['thin_fraction']
            self.options['wmax'] = self.options['emin']*self.options['efrcthn']
            if self.options['wmax'] < 1.0: self.options['wmax'] = 1.0
            self.options['rmax'] = self.options['thin_rmax']
            self.options['thinrat'] = (self.options['emin']*self.options['efrcthn'])/ecut
            self.options['weitrat'] = (self.options['wmax']) # should not be needed, actually
        else:
            logging.critical('Specified thinning method not supported')


    def append_to_steering(self):
        if 'thin' in self.options['other']:
            steering_str  = "THIN   %(efrcthn).2E %(wmax).0f %(rmax).2f      EM thinning level weightmax rmax\n"
            steering_str += "THINH  %(thinrat).2E %(weitrat)f              Ratios for Hadronic thinning\n"
            return steering_str%(self.options)
        return ''

