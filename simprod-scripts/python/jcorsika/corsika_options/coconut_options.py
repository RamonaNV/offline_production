import os, sys

import collections
CoconutOption = collections.namedtuple('CoconutOption', ['name', 'group', 'const', 'help', 'requires', 'conflicts'])
CoconutOption.__new__.__defaults__ = ([], [])
del collections

options_ = {}
options_['75700'] =(
    CoconutOption(*['float_32b',      'float',    '1',  'Force 32bit mode']),
    CoconutOption(*['float_default',  'float',    '2',  'Use compiler default (-m64 on a 64bit machine)']),
    CoconutOption(*['epos',           'he_model', '2',  'EPOS LHC']),
    CoconutOption(*['nexus',          'he_model', '3',  'NEXUS 3.97']),
    CoconutOption(*['qgsjet_i',       'he_model', '4',  'QGSJET 01C']),
    CoconutOption(*['qgsjet_ii',      'he_model', '5',  'QGSJETII-04']),
    CoconutOption(*['sibyll',         'he_model', '6',  'SIBYLL 2.3']),
    CoconutOption(*['venus',          'he_model', '7',  'VENUS 4.12']),
    CoconutOption(*['gheisha',        'le_model', '1',  'GHEISHA 2002d']),
    CoconutOption(*['fluka',          'le_model', '2',  'FLUKA']),
    CoconutOption(*['urqmd',          'le_model', '3',  'URQMD 1.3cr']),
    CoconutOption(*['flat',           'geometry', '1',  'horizontal flat detector array']),
    CoconutOption(*['spherical',      'geometry', '2',  'non-flat (volume) detector geometry']),
    CoconutOption(*['cylindrical',    'geometry', '3',  'vertical string detector geometry']),
    CoconutOption(*['cher_rectangle', 'other',    '1a', 'Cherenkov version for rectangular detector grid']),
    CoconutOption(*['cher_iact',      'other',    '1b', 'Cherenkov version for telescope system (using bernlohr IACT C-routines',
             ['cher_rectangle']]),
    CoconutOption(*['cher_absorption','other',    '1c', 'apply atm. absorption, mirror reflectivity & quantum eff.',
             ['cher_rectangle']]),
    CoconutOption(*['cher_auger',     'other',    '1d', 'Auger Cherenkov longitudinal distribution']),
    CoconutOption(*['trajectory',     'other',    '1e', 'TRAJECTory version to follow motion of source on the sky']),
    CoconutOption(*['lpm',            'other',    '2',  'LPM-effect without thinning']),
    CoconutOption(*['thin',           'other',    '2a', 'THINning version (includes LPM)']),
    CoconutOption(*['multithin',      'other',    '2b', 'Multi-THINning version']),
    CoconutOption(*['preshower',      'other',    '3',  'PRESHOWER version for EeV gammas']),
    CoconutOption(*['neutrino',       'other',    '4',  'NEUTRINO version']),
    CoconutOption(*['nuprim_herwig',  'other',    '4a', 'NUPRIM primary neutrino version with HERWIG']),
    CoconutOption(*['icecube_fifo',   'other',    '4b', 'ICECUBE1 FIFO version']),
    CoconutOption(*['icecube_output', 'other',    '4c', 'ICECUBE2 gzip/pipe output']),
    CoconutOption(*['stack_in',       'other',    '5',  'STACK INput of secondaries, no primary particle']),
    CoconutOption(*['charm',          'other',    '6',  'CHARMed particle/tau lepton version with PYTHIA']),
    CoconutOption(*['tau',            'other',    '6a', 'TAU LEPton version with PYTHIA']),
    CoconutOption(*['slant',          'other',    '7',  'SLANT depth instead of vertical depth for longi-distribution']),
    CoconutOption(*['curved',         'other',    '7a', 'CURVED atmosphere version']),
    CoconutOption(*['upward',         'other',    '7b', 'UPWARD particles version']),
    CoconutOption(*['view_cone',      'other',    '7c', 'view-cone version']),
    CoconutOption(*['shower_plot',    'other',    '8a', 'shower PLOT version (PLOTSH) (only for single events)']),
    CoconutOption(*['shower_plot_2',  'other',    '8b', 'shower PLOT(C) version (PLOTSH2) (only for single events)']),
    CoconutOption(*['ana_hist',       'other',    '8c', 'ANAlysis HISTos & THIN (instead of particle file', ['thin']]),
    CoconutOption(*['auger_hist',     'other',    '8d', 'Auger-histo file & THIN', ['thin', 'auger_info']]),
    CoconutOption(*['muon_hist',      'other',    '8e', 'MUON-histo file']),
    CoconutOption(*['atm_ext',        'other',    '9',  'external atmosphere functions (table interpolation using bernlohr C-routines']),
    CoconutOption(*['atm_efield',     'other',    '9a', 'EFIELD version for electrical field in atmosphere']),
    CoconutOption(*['rigidity',       'other',    '9b', 'RIGIDITY Ooty version rejecting low-energy primaries entering Earth-magnetic field']),
    CoconutOption(*['dyn_stack',      'other',   '10a', 'DYNamic intermediate particle STACK']),
    CoconutOption(*['remote-control', 'other',   '10b', 'Remote Control for Corsika']),
    CoconutOption(*['conex',          'other',    'a',  'CONEX for high energy MC and cascade equations',
             ['thin', 'slant', 'curved', 'upward']]),
    CoconutOption(*['parallel',       'other',    'b',  'PARALLEL treatment of subshowers']),
    CoconutOption(*['coreas',         'other',    'c',  'CoREAS Radio Simulations', ['slant']]),
    CoconutOption(*['inclined',       'other',    'd1', 'Inclined observation plane']),
    CoconutOption(*['root_out',       'other',    'd2', 'ROOT particle OUTput file']),
    CoconutOption(*['coast',          'other',    'd3', 'Use an external COAST user library (COrsika data AccesS Tool',
             ['slant']]),
    CoconutOption(*['int_test',       'other',    'e',  'interaction test version (only for 1st interaction']),
    CoconutOption(*['auger_info',     'other',    'f',  'Auger-info file instead of dbase file']),
    CoconutOption(*['compact',        'other',    'g',  'COMPACT particle output file']),
    CoconutOption(*['muprod',         'other',    'h',  'MUPROD to write decaying muons']),
    CoconutOption(*['history',        'other',    'h2', 'prEHISTORY of muons: mother and grandmother']),
    CoconutOption(*['annitest',       'other',    'k',  'annitest cross-section version (obsolete)' ]),
    CoconutOption(*['hit_auger',      'other',    'l',  'hit Auger detector (steered by AUGSCT)'])
)

options_['75600'] =(
    CoconutOption(*['float_32b',      'float',    '1',  'Force 32bit mode']),
    CoconutOption(*['float_default',  'float',    '2',  'Use compiler default (-m64 on a 64bit machine)']),
    CoconutOption(*['dpmjet',         'he_model', '1',  'DPMJET 2.55']),
    CoconutOption(*['epos',           'he_model', '2',  'EPOS LHC']),
    CoconutOption(*['nexus',          'he_model', '3',  'NEXUS 3.97']),
    CoconutOption(*['qgsjet_i',       'he_model', '4',  'QGSJET 01C']),
    CoconutOption(*['qgsjet_ii',      'he_model', '5',  'QGSJETII-04']),
    CoconutOption(*['sibyll',         'he_model', '6',  'SIBYLL 2.3']),
    CoconutOption(*['venus',          'he_model', '7',  'VENUS 4.12']),
    CoconutOption(*['gheisha',        'le_model', '1',  'GHEISHA 2002d']),
    CoconutOption(*['fluka',          'le_model', '2',  'FLUKA']),
    CoconutOption(*['urqmd',          'le_model', '3',  'URQMD 1.3cr']),
    CoconutOption(*['flat',           'geometry', '1',  'horizontal flat detector array']),
    CoconutOption(*['spherical',      'geometry', '2',  'non-flat (volume) detector geometry']),
    CoconutOption(*['cylindrical',    'geometry', '3',  'vertical string detector geometry']),
    CoconutOption(*['cher_rectangle', 'other',    '1a', 'Cherenkov version for rectangular detector grid']),
    CoconutOption(*['cher_iact',      'other',    '1b', 'Cherenkov version for telescope system (using bernlohr IACT C-routines',
             ['cher_rectangle']]),
    CoconutOption(*['cher_absorption','other',    '1c', 'apply atm. absorption, mirror reflectivity & quantum eff.',
             ['cher_rectangle']]),
    CoconutOption(*['cher_auger',     'other',    '1d', 'Auger Cherenkov longitudinal distribution']),
    CoconutOption(*['trajectory',     'other',    '1e', 'TRAJECTory version to follow motion of source on the sky']),
    CoconutOption(*['lpm',            'other',    '2',  'LPM-effect without thinning']),
    CoconutOption(*['thin',           'other',    '2a', 'THINning version (includes LPM)']),
    CoconutOption(*['multithin',      'other',    '2b', 'Multi-THINning version']),
    CoconutOption(*['preshower',      'other',    '3',  'PRESHOWER version for EeV gammas']),
    CoconutOption(*['neutrino',       'other',    '4',  'NEUTRINO version']),
    CoconutOption(*['nuprim_herwig',  'other',    '4a', 'NUPRIM primary neutrino version with HERWIG']),
    CoconutOption(*['icecube_fifo',   'other',    '4b', 'ICECUBE1 FIFO version']),
    CoconutOption(*['icecube_output', 'other',    '4c', 'ICECUBE2 gzip/pipe output']),
    CoconutOption(*['stack_in',       'other',    '5',  'STACK INput of secondaries, no primary particle']),
    CoconutOption(*['charm',          'other',    '6',  'CHARMed particle/tau lepton version with PYTHIA']),
    CoconutOption(*['tau',            'other',    '6a', 'TAU LEPton version with PYTHIA']),
    CoconutOption(*['slant',          'other',    '7',  'SLANT depth instead of vertical depth for longi-distribution']),
    CoconutOption(*['curved',         'other',    '7a', 'CURVED atmosphere version']),
    CoconutOption(*['upward',         'other',    '7b', 'UPWARD particles version']),
    CoconutOption(*['view_cone',      'other',    '7c', 'view-cone version']),
    CoconutOption(*['shower_plot',    'other',    '8a', 'shower PLOT version (PLOTSH) (only for single events)']),
    CoconutOption(*['shower_plot_2',  'other',    '8b', 'shower PLOT(C) version (PLOTSH2) (only for single events)']),
    CoconutOption(*['ana_hist',       'other',    '8c', 'ANAlysis HISTos & THIN (instead of particle file', ['thin']]),
    CoconutOption(*['auger_hist',     'other',    '8d', 'Auger-histo file & THIN', ['thin', 'auger_info']]),
    CoconutOption(*['muon_hist',      'other',    '8e', 'MUON-histo file']),
    CoconutOption(*['atm_ext',        'other',    '9',  'external atmosphere functions (table interpolation using bernlohr C-routines']),
    CoconutOption(*['atm_efield',     'other',    '9a', 'EFIELD version for electrical field in atmosphere']),
    CoconutOption(*['conex',          'other',    'a',  'CONEX for high energy MC and cascade equations',
             ['thin', 'slant', 'curved', 'upward']]),
    CoconutOption(*['parallel',       'other',    'b',  'PARALLEL treatment of subshowers']),
    CoconutOption(*['coreas',         'other',    'c',  'CoREAS Radio Simulations', ['slant']]),
    CoconutOption(*['inclined',       'other',    'd1', 'Inclined observation plane']),
    CoconutOption(*['root_out',       'other',    'd2', 'ROOT particle OUTput file']),
    CoconutOption(*['coast',          'other',    'd3', 'Use an external COAST user library (COrsika data AccesS Tool',
             ['slant']]),
    CoconutOption(*['int_test',       'other',    'e',  'interaction test version (only for 1st interaction']),
    CoconutOption(*['auger_info',     'other',    'f',  'Auger-info file instead of dbase file']),
    CoconutOption(*['compact',        'other',    'g',  'COMPACT particle output file']),
    CoconutOption(*['muprod',         'other',    'h',  'MUPROD to write decaying muons']),
    CoconutOption(*['history',        'other',    'h2', 'prEHISTORY of muons: mother and grandmother']),
    CoconutOption(*['annitest',       'other',    'k',  'annitest cross-section version (obsolete)' ]),
    CoconutOption(*['hit_auger',      'other',    'l',  'hit Auger detector (steered by AUGSCT)'])
)

options_['75000'] =(
    CoconutOption(*['float_32b',      'float',    '1',  'Force 32bit mode']),
    CoconutOption(*['float_default',  'float',    '2',  'Use compiler default (-m64 on a 64bit machine)']),
    CoconutOption(*['dpmjet',         'he_model', '1',  'DPMJET 2.55']),
    CoconutOption(*['epos',           'he_model', '2',  'EPOS LHC']),
    CoconutOption(*['nexus',          'he_model', '3',  'NEXUS 3.97']),
    CoconutOption(*['qgsjet_i',       'he_model', '4',  'QGSJET 01C']),
    CoconutOption(*['qgsjet_ii',      'he_model', '5',  'QGSJETII-04']),
    CoconutOption(*['sibyll',         'he_model', '6',  'SIBYLL 2.3']),
    CoconutOption(*['venus',          'he_model', '7',  'VENUS 4.12']),
    CoconutOption(*['gheisha',        'le_model', '1',  'GHEISHA 2002d']),
    CoconutOption(*['fluka',          'le_model', '2',  'FLUKA']),
    CoconutOption(*['urqmd',          'le_model', '3',  'URQMD 1.3cr']),
    CoconutOption(*['flat',           'geometry', '1',  'horizontal flat detector array']),
    CoconutOption(*['spherical',      'geometry', '2',  'non-flat (volume) detector geometry']),
    CoconutOption(*['cylindrical',    'geometry', '3',  'vertical string detector geometry']),
    CoconutOption(*['cher_rectangle', 'other',    '1a', 'Cherenkov version for rectangular detector grid']),
    CoconutOption(*['cher_iact',      'other',    '1b', 'Cherenkov version for telescope system (using bernlohr IACT C-routines',
             ['cher_rectangle']]),
    CoconutOption(*['cher_absorption','other',    '1c', 'apply atm. absorption, mirror reflectivity & quantum eff.',
             ['cher_rectangle']]),
    CoconutOption(*['cher_auger',     'other',    '1d', 'Auger Cherenkov longitudinal distribution']),
    CoconutOption(*['trajectory',     'other',    '1e', 'TRAJECTory version to follow motion of source on the sky']),
    CoconutOption(*['thin',           'other',    '2a', 'THINning version (includes LPM)']),
    CoconutOption(*['lpm',            'other',    '2',  'LPM-effect without thinning']),
    CoconutOption(*['preshower',      'other',    '3',  'PRESHOWER version for EeV gammas']),
    CoconutOption(*['neutrino',       'other',    '4',  'NEUTRINO version']),
    CoconutOption(*['nuprim_herwig',  'other',    '4a', 'NUPRIM primary neutrino version with HERWIG']),
    CoconutOption(*['icecube_fifo',   'other',    '4b', 'ICECUBE1 FIFO version']),
    CoconutOption(*['icecube_output', 'other',    '4c', 'ICECUBE2 gzip/pipe output']),
    CoconutOption(*['stack_in',       'other',    '5',  'STACK INput of secondaries, no primary particle']),
    CoconutOption(*['charm',          'other',    '6',  'CHARMed particle/tau lepton version with PYTHIA']),
    CoconutOption(*['tau',            'other',    '6a', 'TAU LEPton version with PYTHIA']),
    CoconutOption(*['slant',          'other',    '7',  'SLANT depth instead of vertical depth for longi-distribution']),
    CoconutOption(*['curved',         'other',    '7a', 'CURVED atmosphere version']),
    CoconutOption(*['upward',         'other',    '7b', 'UPWARD particles version']),
    CoconutOption(*['view_cone',      'other',    '7c', 'view-cone version']),
    CoconutOption(*['shower_plot',    'other',    '8a', 'shower PLOT version (PLOTSH) (only for single events)']),
    CoconutOption(*['shower_plot_2',  'other',    '8b', 'shower PLOT(C) version (PLOTSH2) (only for single events)']),
    CoconutOption(*['ana_hist',       'other',    '8c', 'ANAlysis HISTos & THIN (instead of particle file', ['thin']]),
    CoconutOption(*['auger_hist',     'other',    '8d', 'Auger-histo file & THIN', ['thin', 'auger_info']]),
    CoconutOption(*['muon_hist',      'other',    '8e', 'MUON-histo file']),
    CoconutOption(*['atm_ext',        'other',    '9',  'external atmosphere functions (table interpolation using bernlohr C-routines']),
    CoconutOption(*['atm_efield',     'other',    '9a', 'EFIELD version for electrical field in atmosphere']),
    CoconutOption(*['conex',          'other',    'a',  'CONEX for high energy MC and cascade equations',
             ['thin', 'slant', 'curved', 'upward']]),
    CoconutOption(*['parallel',       'other',    'b',  'PARALLEL treatment of subshowers']),
    CoconutOption(*['coreas',         'other',    'c',  'CoREAS Radio Simulations', ['slant']]),
    CoconutOption(*['inclined',       'other',    'd1', 'Inclined observation plane']),
    CoconutOption(*['root_out',       'other',    'd2', 'ROOT particle OUTput file']),
    CoconutOption(*['coast',          'other',    'd3', 'Use an external COAST user library (COrsika data AccesS Tool',
             ['slant']]),
    CoconutOption(*['int_test',       'other',    'e',  'interaction test version (only for 1st interaction']),
    CoconutOption(*['auger_info',     'other',    'f',  'Auger-info file instead of dbase file']),
    CoconutOption(*['compact',        'other',    'g',  'COMPACT particle output file']),
    CoconutOption(*['muprod',         'other',    'h',  'MUPROD to write decaying muons']),
    CoconutOption(*['history',        'other',    'h2', 'prEHISTORY of muons: mother and grandmother']),
    CoconutOption(*['annitest',       'other',    'k',  'annitest cross-section version (obsolete)' ]),
    CoconutOption(*['hit_auger',      'other',    'l',  'hit Auger detector (steered by AUGSCT)'])
)

options_['74005'] =(
    CoconutOption(*['float_32b',      'float',    '1',  'Force 32bit mode']),
    CoconutOption(*['float_default',  'float',    '2',  'Use compiler default (-m64 on a 64bit machine)']),
    CoconutOption(*['dpmjet',         'he_model', '1',  'DPMJET 2.55']),
    CoconutOption(*['epos',           'he_model', '2',  'EPOS LHC']),
    CoconutOption(*['nexus',          'he_model', '3',  'NEXUS 3.97']),
    CoconutOption(*['qgsjet_i',       'he_model', '4',  'QGSJET 01C']),
    CoconutOption(*['qgsjet_ii',      'he_model', '5',  'QGSJETII-04']),
    CoconutOption(*['sibyll',         'he_model', '6',  'SIBYLL 2.1']),
    CoconutOption(*['venus',          'he_model', '7',  'VENUS 4.12']),
    CoconutOption(*['gheisha',        'le_model', '1',  'GHEISHA 2002d']),
    CoconutOption(*['fluka',          'le_model', '2',  'FLUKA']),
    CoconutOption(*['urqmd',          'le_model', '3',  'URQMD 1.3cr']),
    CoconutOption(*['flat',           'geometry', '1',  'horizontal flat detector array']),
    CoconutOption(*['spherical',      'geometry', '2',  'non-flat (volume) detector geometry']),
    CoconutOption(*['cylindrical',    'geometry', '3',  'vertical string detector geometry']),
    CoconutOption(*['cher_rectangle', 'other',    '1a', 'Cherenkov version for rectangular detector grid']),
    CoconutOption(*['cher_iact',      'other',    '1b', 'Cherenkov version for telescope system (using bernlohr IACT C-routines',
             ['cher_rectangle']]),
    CoconutOption(*['cher_absorption','other',    '1c', 'apply atm. absorption, mirror reflectivity & quantum eff.',
             ['cher_rectangle']]),
    CoconutOption(*['cher_auger',     'other',    '1d', 'Auger Cherenkov longitudinal distribution']),
    CoconutOption(*['trajectory',     'other',    '1e', 'TRAJECTory version to follow motion of source on the sky']),
    CoconutOption(*['thin',           'other',    '2',  'THINning version']),
    CoconutOption(*['lpm',            'other',    '2e', 'LPM-effect without thinning']),
    CoconutOption(*['preshower',      'other',    '3',  'PRESHOWER version for EeV gammas']),
    CoconutOption(*['neutrino',       'other',    '4',  'NEUTRINO version']),
    CoconutOption(*['nuprim_herwig',  'other',    '4a', 'NUPRIM primary neutrino version with HERWIG']),
    CoconutOption(*['stack_in',       'other',    '5',  'STACK INput of secondaries, no primary particle']),
    CoconutOption(*['charm',          'other',    '6',  'CHARMed particle/tau lepton version with PYTHIA']),
    CoconutOption(*['tau',            'other',    '6a', 'TAU LEPton version with PYTHIA']),
    CoconutOption(*['slant',          'other',    '7',  'SLANT depth instead of vertical depth for longi-distribution']),
    CoconutOption(*['curved',         'other',    '7a', 'CURVED atmosphere version']),
    CoconutOption(*['upward',         'other',    '7b', 'UPWARD particles version']),
    CoconutOption(*['view_cone',      'other',    '7c', 'view-cone version']),
    CoconutOption(*['shower_plot',    'other',    '8a', 'shower PLOT version (PLOTSH (only for single events']),
    CoconutOption(*['shower_plot_2',  'other',    '8b', 'shower PLOT(C version (PLOTSH2 (only for single events']),
    CoconutOption(*['ana_hist',       'other',    '8c', 'ANAlysis HISTos & THIN (instead of particle file',['thin']]),
    CoconutOption(*['auger_hist',     'other',    '8d', 'Auger-histo file & THIN', ['thin', 'auger_info']]),
    CoconutOption(*['muon_hist',      'other',    '8e', 'MUON-histo file']),
    CoconutOption(*['atm_ext',        'other',    '9',  'external atmosphere functions (table interpolation using bernlohr C-routines']),
    CoconutOption(*['atm_efield',     'other',    '9a', 'EFIELD version for electrical field in atmosphere']),
    CoconutOption(*['conex',          'other',    'a',  'CONEX for high energy MC and cascade equations',
             ['thin', 'slant', 'curved', 'upward']]),
    CoconutOption(*['parallel',       'other',    'b',  'PARALLEL treatment of subshowers']),
    CoconutOption(*['coreas',         'other',    'c',  'CoREAS Radio Simulations', ['slant']]),
    CoconutOption(*['inclined',       'other',    'd1', 'Inclined observation plane']),
    CoconutOption(*['root_out',       'other',    'd2', 'ROOT particle OUTput file', ['slant']]),
    CoconutOption(*['coast',          'other',    'd3', 'Use an external COAST user library (COrsika data AccesS Tool']),
    CoconutOption(*['int_test',       'other',    'e',  'interaction test version (only for 1st interaction']),
    CoconutOption(*['auger_info',     'other',    'f',  'Auger-info file instead of dbase file']),
    CoconutOption(*['compact',        'other',    'g',  'COMPACT particle output file']),
    CoconutOption(*['muprod',         'other',    'h',  'MUPROD to write decaying muons']),
    CoconutOption(*['history',        'other',    'h2', 'prEHISTORY of muons: mother and grandmother'])
)


options_['73700'] =(
    CoconutOption(*['float_32b',      'float',    '1',  'Force 32bit mode']),
    CoconutOption(*['float_default',  'float',    '2',  'Use compiler default (-m64 on a 64bit machine)']),
    CoconutOption(*['dpmjet',         'he_model', '1',  'DPMJET 2.55']),
    CoconutOption(*['epos',           'he_model', '2',  'EPOS LHC']),
    CoconutOption(*['nexus',          'he_model', '3',  'NEXUS 3.97']),
    CoconutOption(*['qgsjet_i',       'he_model', '4',  'QGSJET 01C']),
    CoconutOption(*['qgsjet_ii',      'he_model', '5',  'QGSJETII-04']),
    CoconutOption(*['sibyll',         'he_model', '6',  'SIBYLL 2.1']),
    CoconutOption(*['venus',          'he_model', '7',  'VENUS 4.12']),
    CoconutOption(*['gheisha',        'le_model', '1',  'GHEISHA 2002d']),
    CoconutOption(*['fluka',          'le_model', '2',  'FLUKA']),
    CoconutOption(*['urqmd',          'le_model', '3',  'URQMD 1.3cr']),
    CoconutOption(*['flat',           'geometry', '1',  'horizontal flat detector array']),
    CoconutOption(*['spherical',      'geometry', '2',  'non-flat (volume) detector geometry']),
    CoconutOption(*['cylindrical',    'geometry', '3',  'vertical string detector geometry']),
    CoconutOption(*['cher_rectangle', 'other',    '1', 'Cherenkov version for rectangular detector grid']),
    CoconutOption(*['cher_iact',      'other',    '2', 'Cherenkov version for telescope system (using bernlohr IACT C-routines',
             ['cher_rectangle']]),
    CoconutOption(*['cher_absorption','other',    '3', 'apply atm. absorption, mirror reflectivity & quantum eff.',
             ['cher_rectangle']]),
    CoconutOption(*['atm_ext',        'other',    '4',  'external atmosphere functions (table interpolation using bernlohr C-routines']),
    CoconutOption(*['cher_auger',     'other',    'g', 'Auger Cherenkov longitudinal distribution']),
    CoconutOption(*['trajectory',     'other',    'u', 'TRAJECTory version to follow motion of source on the sky']),
    CoconutOption(*['thin',           'other',    '5',  'THINning version']),
    CoconutOption(*['lpm',            'other',    'l', 'LPM-effect without thinning']),
    CoconutOption(*['preshower',      'other',    'h',  'PRESHOWER version for EeV gammas']),
    CoconutOption(*['neutrino',       'other',    '6',  'NEUTRINO version']),
    CoconutOption(*['nuprim_herwig',  'other',    'n', 'NUPRIM primary neutrino version with HERWIG']),
    CoconutOption(*['stack_in',       'other',    'm',  'STACK INput of secondaries, no primary particle']),
    CoconutOption(*['charm',          'other',    'q',  'CHARMed particle/tau lepton version with PYTHIA']),
    CoconutOption(*['tau',            'other',    'qt', 'TAU LEPton version with PYTHIA']),
    CoconutOption(*['slant',          'other',    '9',  'SLANT depth instead of vertical depth for longi-distribution']),
    CoconutOption(*['curved',         'other',    'a', 'CURVED atmosphere version']),
    CoconutOption(*['upward',         'other',    'b', 'UPWARD particles version']),
    CoconutOption(*['view_cone',      'other',    'c', 'view-cone version']),
    CoconutOption(*['shower_plot',    'other',    '7', 'shower PLOT version (PLOTSH (only for single events']),
    CoconutOption(*['shower_plot_2',  'other',    '72', 'shower PLOT(C version (PLOTSH2 (only for single events']),
    CoconutOption(*['ana_hist',       'other',    'd', 'ANAlysis HISTos & THIN (instead of particle file', ['thin']]),
    CoconutOption(*['auger_hist',     'other',    'f', 'Auger-histo file & THIN', ['thin', 'auger_info']]),
    CoconutOption(*['atm_efield',     'other',    'v', 'EFIELD version for electrical field in atmosphere']),
    CoconutOption(*['conex',          'other',    'w',  'CONEX for high energy MC and cascade equations',
             ['thin', 'slant', 'curved', 'upward']]),
    CoconutOption(*['parallel',       'other',    'p',  'PARALLEL treatment of subshowers']),
    CoconutOption(*['root_out',       'other',    'o', 'ROOT particle OUTput file']),
    CoconutOption(*['coast',          'other',    't', 'Use an external COAST user library (COrsika data AccesS Tool',
             ['root_out', 'slant']]),
    CoconutOption(*['int_test',       'other',    '8',  'interaction test version (only for 1st interaction']),
    CoconutOption(*['auger_info',     'other',    'e',  'Auger-info file instead of dbase file']),
    CoconutOption(*['compact',        'other',    'j',  'COMPACT particle output file']),
    CoconutOption(*['muprod',         'other',    'i',  'MUPROD to write decaying muons']),
    CoconutOption(*['history',        'other',    's', 'prEHISTORY of muons: mother and grandmother'])
)


class CoconutOptionCollection(object):
    """
    Class to handle the conversion from a human-readable option format
    (for optparse) and the options to be used with coconut. The actual options
    in coconut are version-dependent, so one needs to create an instance of this
    class for each version.

    The available versions are stored in the available_versions data member.

    
    """
    options = None
    name_from_option = None
    option_from_name = None

    available_versions = sorted(options_.keys())

    def __init__(self, version):
        """
        Initialize the options according to a given CORSIKA version.
        If no version is given, it will try to get it from the command
        line option --version.
        """
        version = str(version)
        if not version:
            raise Exception('CORSIKA version was not specified')
        if not version in options_.keys():
            raise Exception('%s is not one of the available CORSIKA versions: %s'%(version, ', '.join(options_.keys())))
        self.options = options_[version]
        self.name_from_option = {(o.group, o.const):o.name for o in self.options}
        self.option_from_name = {c:(a,b) for (a,b),c in self.name_from_option.iteritems()}
        self.option_groups = []
        for g in [o.group for o in self.options]:
            if not g in self.option_groups: self.option_groups.append(g)
        self.options = dict([(o.name, o) for o in options_[version]])

    def get_coconut_options(self, options):
        """
        Convert a dictionary of options (in a readable format) to the actual option values used in coconut.
        The options in coconut change with versions.

        This function converts something like this (from option parsing or whatever):
            { 'float': 'float_default',
              'he_model': 'sibyll',
              'le_model': 'fluka',
              'geometry': 'flat',
              'other': ['slant', 'neutrino', 'muprod', 'history']
            }
        into something like this (actual coconut options):
            { 'float': '2',
              'geometry': '1',
              'he_model': '6',
              'le_model': '2',
              'other': ['9', '6', 'i', 's']
            }
        Note that each key corresponds to a section in coconut.
        """
        for n,v in options.iteritems():
            if n in ['float', 'he_model', 'le_model', 'geometry'] and not v in self.option_from_name.keys():
                raise Exception('%s "%s" is not an available option'%(n,v))
        opt = dict([self.option_from_name[v] for n,v in options.iteritems()
                    if n in ['float', 'he_model', 'le_model', 'geometry']])
        opt['other'] = [self.option_from_name[n][1] for n in options['other']] if 'other' in options.keys() else []
        # checking that there is one option for each group and one group for each option
        assert(sorted(opt.keys()) == sorted(set([k[0] for k in self.name_from_option.keys()])))
        return opt

    def name_from_options(self, opt, names=False):
        """
        Name of an executable from the compile options.

        Input options like this:
            { 'float': '2',
              'geometry': '1',
              'he_model': '6',
              'le_model': '2',
              'other': ['9', '6', 'i', 's']
            }
        result in names like this:
            '2_6_2_1_6_9_i_s'

        if names==True, the output looks more like this:
            'float_default_sibyll_fluka_flat_neutrino_slant_muprod_history'
        """
        if 'other' in opt and opt['other']:
            other = sorted(opt['other'], key=lambda o: self.options[self.name_from_option[('other',o)]].const)
            for o in other:
                for r in self.options[self.name_from_option[('other',o)]].requires:
                    if r in other: other.remove(r)
            other = [('other', self.options[self.name_from_option[('other',o)]]) for o in other]
        else:
            other = []
        opt = [(k,self.options[self.name_from_option[(k,opt[k])]]) for k in ['float', 'he_model', 'le_model', 'geometry']] + other
        if names:
            return '_'.join([n[1].name for n in opt])
        return '_'.join([n[1].const for n in opt])

    def options_from_name(self, name, names=False):
        """
        Compile options of an executable from its name.

        Input like this:
            '2_6_2_1_6_9_i_s'
        result in options like this:
            { 'float': '2',
              'geometry': '1',
              'he_model': '6',
              'le_model': '2',
              'other': ['9', '6', 'i', 's']
            }

        if names==True, the output looks more like this:
            { 'float': 'float_default',
              'geometry': 'flat',
              'he_model': 'sibyll',
              'le_model': 'fluka',
              'other': ['neutrino', 'slant', 'muprod', 'history']
            }
        """
        fields = name.split('_')
        if names:
            opt = {k:self.name_from_option[(k,v)] for k,v in zip(['float', 'he_model', 'le_model', 'geometry'], fields[:4])}
            opt['other'] = [self.name_from_option[('other',v)] for v in fields[4:]]
        else:
            opt = dict(zip(['float', 'he_model', 'le_model', 'geometry'], fields[:4]))
            opt['other'] = fields[4:]
        return opt

    def old_name_from_options(self, options):
        """
        Name of an executable from the compile options using the old convention.

        Older tarballs have executables with names like these:
          corsika73700Linux_EPOS_fluka_curved
          corsika73700Linux_EPOS_fluka_thin
          corsika73700Linux_EPOS_fluka_curved_thin
          corsika73700Linux_EPOS_fluka

        Note that the names of the models do not correspond to my "standard" names
        """
        names = { 'dpmjet'    :'DPMJET',
                  'epos'      :'EPOS',
                  'nexus'     :'NEXUS',
                  'qgsjet_i'  :'QGSJET',
                  'qgsjet_ii' :'QGSII',
                  'sibyll'    :'SIBYLL',
                  'venus'     :'VENUS',
                  'gheisha'   :'gheisha',
                  'fluka'     :'fluka',
                  'urqmd'     :'urqmd'
        }

        name = names[self.name_from_option[('he_model',options['he_model'])]]
        name += '_' + names[self.name_from_option[('le_model',options['le_model'])]]
        if self.option_from_name['curved'][1] in options['other']:
            name += '_curved'
        if self.option_from_name['thin'][1] in options['other']:
            name += '_thin'
        return name
