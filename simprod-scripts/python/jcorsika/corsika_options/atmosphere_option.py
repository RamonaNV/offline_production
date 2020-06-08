"""
 copyright  (c) 2005 the icecube collaboration

 @version: $Revision: $
 @date: $Date: $
 @author: Javier Gonzalez <javierg@udel.edu>
 """

from option import OptionSuite, Option, OptionGroup
import logging

class AtmosphereOption(OptionSuite):
    """
    Option to set from a list of atmospheric profiles that are not in CORSIKA.
    """
    def __init__(self, options, groups):
        OptionSuite.__init__(self, options)

        options.append(Option(name='atmod',
                              long_name='--atmod',
                              help = 'Atmosphere model (October=13)',
                              default = 13,
                              type =  'int'))

        options.append(Option(name='ratmo',
                              long_name='--ratmo',
                              help = 'Use a "real" atmosphere. If this is set, atmod is ignored.',
                              type =  'int'))


    def append_to_steering(self):
        ratmo = self.options['ratmo']
        steering_str = ''
        if ratmo == 1:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -91.6956   7.01491   0.505452  -0.00181302   0.00207722\n"
            steering_str +="ATMB        1125.71    1149.81   1032.68   490.789\n"
            steering_str +="ATMC        821621     635444    682968    807327        5.4303203E9\n"
            steering_str +="ATMLAY      780000     1640000   4040000   10000000\n"
        elif ratmo == 2:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -72.1988   22.7002   0.430171  -0.001203     0.00207722\n"
            steering_str +="ATMB        1108.19    1159.77   1079.25   523.956\n"
            steering_str +="ATMC        786271     599986    667432    780919        5.4303203E9\n"
            steering_str +="ATMLAY      800000     1060000   4040000   10000000\n"
        elif ratmo == 3:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -63.729    -1.02799  0.324414  -0.000490772  0.00207722\n"
            steering_str +="ATMB        1102.66    1093.56   1198.93   589.827\n"
            steering_str +="ATMC        764831     660389    636118    734909        5.4303203E9\n"
            steering_str +="ATMLAY      670000     2240000   4040000   10000000\n"
        elif ratmo == 4:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -69.7259   -2.79781  0.262692  -8.41695e-05  0.00207722\n"
            steering_str +="ATMB        1111.7     1128.64   1413.98   587.688\n"
            steering_str +="ATMC        766099     641716    588082    693300        5.4303203E9\n"
            steering_str +="ATMLAY      760000     2200000   4040000   10000000\n"
        elif ratmo == 5:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -78.5551   -5.33239  0.312889  -9.20472e-05  0.00152236\n"
            steering_str +="ATMB        1118.46    1169.09   1577.71   452.177\n"
            steering_str +="ATMC        776648     626683    553087    696835        7.4095699E9\n"
            steering_str +="ATMLAY      840000     2000000   3970000   10000000\n"
        elif ratmo == 6:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -92.6125   -8.5645   0.363986  1.65164e-05   0.00207722\n"
            steering_str +="ATMB        1129.88    1191.98   1619.82   411.586\n"
            steering_str +="ATMC        791177     618840    535235    692253        5.4303203E9\n"
            steering_str +="ATMLAY      850000     1790000   3840000   10000000\n"
        elif ratmo == 7:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA       -89.9639    -13.9697  0.441631  -1.46525e-05  0.00207722\n"
            steering_str +="ATMB        1125.73    1180.47   1581.43   373.796\n"
            steering_str +="ATMC        784553     628042    531652    703417        5.4303203E9\n"
            steering_str +="ATMLAY      850000     1590000   3750000   10000000\n"
        elif ratmo == 8:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -90.4253   -18.7154  0.51393   -0.00021565   0.00152236\n"
            steering_str +="ATMB        1125.01    1175.6    1518.03   299.006\n"
            steering_str +="ATMC        781628     633793    533269    737794        7.4095699E9\n"
            steering_str +="ATMLAY      850000     1440000   3750000   10000000\n"
        elif ratmo == 9:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -91.686    -23.3519  0.891302  -0.000765666  0.00207722\n"
            steering_str +="ATMB        1125.53    1169.77   1431.26   247.03\n"
            steering_str +="ATMC        786017     645241    545022    805419        5.4303203E9\n"
            steering_str +="ATMLAY      850000     1300000   3620000   10000000\n"
        elif ratmo == 10:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        451.616    -85.5456  2.06082   -0.001076     0.00207722\n"
            steering_str +="ATMB        849.239    1113.16   1322.28   372.242\n"
            steering_str +="ATMC        225286     789340    566132    796434        5.4303203E9\n"
            steering_str +="ATMLAY      310000     1010000   3150000   10000000\n"
        elif ratmo == 11:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -152.853   4.22741   1.38352   -0.00115014   0.00207722\n"
            steering_str +="ATMB        1174.09    1272.49   975.906   481.615\n"
            steering_str +="ATMC        891602     582119    643130    783786        5.4303203E9\n"
            steering_str +="ATMLAY      850000     2240000   3240000   10000000\n"
        elif ratmo == 12:
            steering_str += "ATMOD      10                              real atmosphere\n"
            steering_str +="ATMA        -100.386   5.43849   0.399465  -0.00175472   0.00207722\n"
            steering_str +="ATMB        1128.71    1198.1    858.522   480.142\n"
            steering_str +="ATMC        829352     612649    706104    806875        5.4303203E9\n"
            steering_str +="ATMLAY      850000     2200000   4040000   10000000\n"
        elif ratmo == 13:
            steering_str += "ATMOD       10                              real atmosphere\n"
            steering_str += "ATMA        -69.7259   -2.79781  0.262692  -8.41695e-05 0.00207722\n"
            steering_str += "ATMB        1111.7     1128.64   1413.98   587.688\n"
            steering_str += "ATMC        766099     641716    588082    693300       5.4303203E9\n"
            steering_str += "ATMLAY      760000     2200000   4040000   10000000\n"
        else:
            steering_str += "ATMOD   %(atmod)s                         october atmosphere\n"

        return steering_str%(self.options)
