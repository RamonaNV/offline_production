"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
 """

from option import OptionSuite, Option, OptionGroup
import logging

class ParametersOption(OptionSuite):
    """
    Main Parameters that are usually set.
    """
    def __init__(self, options, groups):
        OptionSuite.__init__(self, options)

        groups.append(OptionGroup(name='other',
                                  title='Other Options',
                                  help='Less common options'))

        options.append(Option('runnr',
                              long_name='--runnr',
                              help='Run number',
                              default = 1,
                              type =  'int'))
        
        options.append(Option(long_name='--evtnr',
                              name = 'evtnr',
                              help = 'Number of first event',
                              default = 1,
                              type =  'int'))

        options.append(Option(long_name='--seed',
                              name = 'seed',
                              help = 'Random seed',
                              default = None,
                              type =  'int'))

        options.append(Option(long_name='--egs_seed_offset' ,
                              name = 'egs_seed_offset' ,
                              help = 'value to be added to EGS seed (for debugging)',
                              default = 0,
                              type = 'int',
                              group='other'))

        options.append(Option(long_name='--nshow',
                              name = 'nshow',
                              help = 'Number of Showers',
                              default = 1,
                              type = 'int'))

        options.append(Option(long_name='--no-nkg',
                              name = 'donkg',
                              help = 'Run NKG',
                              default = True,
                              action = "store_false",
                              group='other'))

        options.append(Option(long_name='--no-egs',
                              name = 'doegs',
                              help = 'Run EGS',
                              default = True,
                              action = "store_false",
                              group='other'))

        options.append(Option(long_name='--eslope',
                              name = 'eslope',
                              help = 'CR spectral index (only if ranpri=0)',
                              default = -2.7,
                              type =  'float'))

        options.append(Option(long_name='--prmpar',
                              name = 'prmpar',
                              help = 'CR primary type',
                              default = 14,
                              type =  'int'))

        options.append(Option(long_name='--theta-min',
                              name = 'cthmin',
                              help = 'Min theta of injected cosmic rays',
                              default = 0.0,
                              type =  'float'))

        options.append(Option(long_name='--theta-max',
                              name = 'cthmax',
                              help = 'Max theta of injected cosmic rays',
                              default = 65.,
                              type =  'float'))

        options.append(Option(long_name='--phi-min',
                              name = 'cphmin',
                              help = 'Min phi of injected cosmic rays',
                              default = 0.0,
                              type =  'float'))

        options.append(Option(long_name='--phi-max',
                              name = 'cphmax',
                              help = 'Max phi of injected cosmic rays',
                              default = 360.0,
                              type =  'float'))

        options.append(Option(long_name='--emin',
                              name = 'emin',
                              help = 'CR min energy',
                              default = 600.,
                              type =  'float'))

        options.append(Option(long_name='--emax',
                              name = 'emax',
                              help = 'CR max energy',
                              default = 1.e11,
                              type =  'float'))

        options.append(Option(long_name='--ecuts1',
                              name = 'ecuts1',
                              help = 'hadron min energy (see corsika docs)',
                              default = 0.1,
                              type =  'float',
                              group='other'))

        options.append(Option(long_name='--ecuts2',
                              name = 'ecuts2',
                              help = 'muon min energy (see corsika docs)',
                              default = 0.1,
                              type =  'float',
                              group='other'))

        options.append(Option(long_name='--ecuts3',
                              name = 'ecuts3',
                              help = 'electron min energy (see corsika docs)',
                              default = 0.00025,
                              type =  'float',
                              group='other'))

        options.append(Option(long_name='--ecuts4',
                              name = 'ecuts4',
                              help = 'photon min energy (see corsika docs)',
                              default = 0.00025,
                              type =  'float',
                              group='other'))

        options.append(Option(long_name='--nuaddi',
                              name = 'nuaddi',
                              help = 'additional information for neutrinos',
                              default = False,
                              action = "store_true",
                              group='other'))

        options.append(Option(long_name='--muaddi',
                              name = 'muaddi',
                              help = 'additional information for muons',
                              default = False,
                              action = "store_true",
                              group='other'))

        options.append(Option(long_name='--emaddi',
                              name = 'emaddi',
                              help = 'additional information for EM particles',
                              default = False,
                              action = "store_true",
                              group='other'))

        options.append(Option(long_name='--obslev',
                              name = 'obslev',
                              help = 'distance above sea level (in cm)',
                              type = float,
                              default = 2837.e2))

        options.append(Option(long_name='--length',
                              name = 'length',
                              help = 'length of generation cylinder in m '
                              '(for detcfg = length/2*radius calculation)',
                              default = 1400.,
                              type =  'float'))

        options.append(Option(long_name='--radius',
                              name = 'radius'          ,
                              help = 'radius of generation cylinder in m '
                              '(for detcfg = length/2*radius calculation)',
                              default = 700.,
                              type =  'float'))

        options.append(Option(long_name='--kcut',
                              name = 'kcut',
                              help = 'minimum neutrino energy required to keep the shower',
                              default = None,
                              type =  'float',
                              group='other'))

        options.append(Option(long_name='--save-long',
                              name = 'save_long',
                              help = 'Save the longitudinal profile in output file (LONG blocks).',
                              default = False,
                              action = "store_true",
                              group='other'))

        options.append(Option(long_name='--long-step',
                              name = 'long_step',
                              help = 'Step size in longitudinal distributions',
                              default = 20.,
                              type =  'float',
                              group='other'))

        options.append(Option(long_name='--arrang',
                              name = 'arrang',
                              help = 'Rotation of detector relative to magnetic '
                              'field direction (in degrees)',
                              default = 0.,
                              type = 'float',
                              group='other'))

        options.append(Option(long_name='--debug',
                              name = 'debug',
                              help = 'Enable disable debug mode',
                              default = False,
                              action = "store_true",
                              group='other'))

        options.append(Option(long_name='--hi-low',
                              name = 'hi_low',
                              help = 'Transition Energy between Models',
                              default = None,
                              type =  'float',
                              group='other'))

        options.append(Option(long_name='--mag-vertical',
                              name = 'mag_vertical',
                              help = 'Magnetic field vertical component (micro-Tesla). The default is for South Pole.',
                              default = 16.4,
                              type =  'float'))

        options.append(Option(long_name='--mag-horizontal',
                              name = 'mag_horizontal',
                              help = 'Magnetic field horizontal component (micro-Tesla). The default is for South Pole.',
                              default = -53.4,
                              type =  'float'))

        options.append(Option(long_name='--fix-height',
                              name = 'fixhei',
                              help = 'Fix height of first interaction.',
                              default = None,
                              type =  'float',
                              group='other'))

        options.append(Option(long_name='--fixchi',
                              name = 'fixchi',
                              help = 'Starting altitude of primary particle.',
                              default = None,
                              type =  'float',
                              group='other'))


    def append_to_steering(self):
        bool_str = lambda v: 'T' if v else 'F'

        # Convert input phi from IceCube Coordinates to CORSIKA Coordinates.
        # CORSIKA will rotate the particles back to IceCube Coordinates in
        # the output routine.  Also, IceCube Zenith and Azimuth point to
        # where the particle came from, while CORSIKA points to where it
        # is going.  Also CORSIKA measures zenith from -z.
        cphmin_cc = (self.options['cphmin'] + self.options['arrang'] + 180)
        cphmax_cc = (self.options['cphmax'] + self.options['arrang'] + 180)
        while cphmin_cc > cphmax_cc:
            cphmax_cc += 360
        while cphmax_cc - cphmin_cc > 360:
            cphmax_cc -= 360
        while cphmin_cc > 0 and cphmax_cc > 0:
            cphmin_cc -= 360
            cphmax_cc -= 360

        if not self.options['seed']: self.options['seed'] = 3*self.options['runnr']

        steering_str = ""

        steering_str += "RUNNR %d \n" % self.options['runnr']
        steering_str += "EVTNR %d \n" % self.options['evtnr']
        steering_str += "NSHOW %d \n" % self.options['nshow']
        steering_str += "PRMPAR %d \n" % self.options['prmpar']
        steering_str += "ESLOPE %d \n" % self.options['eslope']
        steering_str += "ERANGE %f %f \n" % (self.options['emin'], self.options['emax'])
        steering_str += "THETAP %f %f \n" % (self.options['cthmin'], self.options['cthmax'])
        steering_str += "PHIP %f %f \n" % (cphmin_cc, cphmax_cc)

        steering_str += "SEED %d 0 0 \n" % self.options['seed']
        steering_str += "SEED %d 0 0 \n" % (self.options['seed'] + 1 + self.options['egs_seed_offset'])
        steering_str += "SEED %d 0 0 \n" % (self.options['seed'] + 2)
        steering_str += "OBSLEV %.0f \n" % self.options['obslev']

        steering_str += "ELMFLG %s %s \n" %( bool_str(self.options['donkg']), bool_str(self.options['doegs']) )
        if self.options['donkg']:
            steering_str += "RADNKG 2.E5 \n"

        steering_str += "MAGNET  {mag1:.2f} {mag2:.2f} \n".format(mag1=self.options['mag_vertical'], mag2=self.options['mag_horizontal'])

        steering_str += "LONGI T %f T %s \n"%(self.options['long_step'], bool_str(self.options['save_long']))

        steering_str += "HADFLG  0  1  0  1  0  2 \n"
        steering_str += "ECUTS  %.04f %.04f %.04f %.04f \n" % \
            (self.options['ecuts1'], self.options['ecuts2'],
             self.options['ecuts3'], self.options['ecuts4'])

        if self.options['arrang']:
            steering_str += "ARRANG %f \n" % self.options['arrang']

        if self.options['fixhei']:
            steering_str += "FIXHEI  %f  0 \n"%self.options['fixhei']
        if self.options['fixchi']:
            steering_str += "FIXCHI  %f \n"%self.options['fixchi']

        if self.options['emaddi']:
            steering_str += "EMADDI T \n"
        if self.options['muaddi']:
            steering_str += "MUADDI T \n"
        if self.options['nuaddi']:
            steering_str += "NUADDI T \n"

        if self.options['kcut']:
            steering_str += "KCUT %s \n" % self.options['kcut']

        if self.options['hi_low']:
            steering_str += "HILOW   %f"%self.options['hi_low']

        steering_str += "MAXPRT 0 \n"   # maximum number of events with detailed printout
        steering_str += "ECTMAP 100 \n" # cut on gamma factor above which EM particles are printed to MONIOU when passing observing levels
        steering_str += "STEPFC 1.0 \n" # this is corsika (multiple scattering step length factor)
        steering_str += "DEBUG %s 6 F 1000000 \n" % bool_str(self.options['debug'])

        return steering_str

    def configure(self, options):
        # call the base class method
        # this is not 'pure virtual' in the parlance of our times
        OptionSuite.configure(self, options)
        if options['geometry'] == 'flat':
            # this conforms to the behavior of the old simprod code
            options['length'] = None
            options['radius'] = None
            options['detcfg'] = None
        else:
            options['detcfg'] = options['length']/(2.*options['radius'])

    def validate(self):
        ok = True
        if (self.options['cthmin'] > 65 or self.options['cthmax'] > 65) and not 'curved' in self.options['other']:
            ok = False
            logging.warning('Trying to simulate between %.2f and %.2f degree zenith angle with a flat atmosphere is not a good idea.'%(self.options['cthmin'], self.options['cthmax']))
        if self.options['emin'] > self.options['emax']:
            ok = False
            logging.warning('Emin (%.2f) is larger than Emax (%.2f).'%(self.options['emin'], self.options['emax']))
        return ok
