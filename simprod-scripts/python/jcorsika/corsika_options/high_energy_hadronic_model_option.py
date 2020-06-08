"""
copyright  (c) 2005 the icecube collaboration

@version: $Revision: $
@date: $Date: $
@author: Alex Olivas <olivas@icecube.umd.edu> and Javier Gonzalez <javierg@udel.edu>
"""

from option import OptionSuite, Option, OptionGroup
import logging

class HighEnergyHadronicModelOption(OptionSuite):
    """
    FIXME : Write docs here.
    """
    def __init__(self, options, groups):
        OptionSuite.__init__(self, options)

        # this is set in the compile options part
        #self.option_parser.add_option('--he_model',
        #                                dest = 'he_model',
        #                                help = 'Hadronic interaction model',
        #                                default = 'SIBYLL',
        #                                action = "store")

    def append_to_steering(self):

        steering_str = ""

        model = self.options['he_model'].lower()
        if model in ("qgsjet","qgsjet_ii"):
            steering_str += "QGSJET  T  0 \n"
            steering_str += "QGSSIG   T \n"
        elif model == "dpmjet":
            steering_str += "DPMJET  T  0 \n"
            steering_str += "DPJSIG  T \n"
        elif model == "sibyll":
            steering_str += "SIBYLL  T  0 \n"
            steering_str += "SIBSIG  T \n"
        elif model == "epos":
            steering_str += "EPOS    T  0 \n"
            steering_str += "EPOSIG  T \n"
            steering_str += "EPOPAR input epos.param \n"
            steering_str += "EPOPAR fname inics epos.inics \n"
            steering_str += "EPOPAR fname iniev epos.iniev \n"
            steering_str += "EPOPAR fname initl epos.initl \n"
            steering_str += "EPOPAR fname inirj epos.inirj \n"
            steering_str += "EPOPAR fname inihy epos.ini1b \n"
            steering_str += "EPOPAR fname check none \n"
            steering_str += "EPOPAR fname histo none \n"
            steering_str += "EPOPAR fname data  none \n"
            steering_str += "EPOPAR fname copy  none \n"
        #elif model == "hdpm":
        #    steering_str += "HADFLG  0  1  0  1  0  2 \n"
        else:
            raise Exception('Unknown HE model %s' % model)


        return steering_str


    def validate(self):
        #also done in compilation_options
        return True
