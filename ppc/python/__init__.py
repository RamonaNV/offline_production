from __future__ import print_function
from icecube import dataclasses

from icecube.load_pybindings import load_pybindings
load_pybindings(__name__,__path__)

#helper function to build the tau_dnde vector needed /fhl
def tau_dnde_builder(it):
    tmp = []
    for i in it:
        if len(i) != 2:
            raise ValueError("All elements must have exactly two parameters (tau and dN/dE)!")
        if i[0] <= 0:
            raise ValueError("tau must be greater than 0")
        tmp.append(dataclasses.make_pair(i[0],i[1]))
    return dataclasses.I3VectorDoubleDouble(tmp)

def MakeCLSimPropagator(DetectorParams, UseCPUs=False, UseGPUs=True):
    """
    Configure ppc to for use with clsim. 
    
    :param DetectorParams: configuration dict returned by clsim.traysegments.common.setupDetector()
    :param UseCPUs: use CPUs
    :param UseGPUs: use GPUs
    """
    from icecube.icetray import I3Units
    from icecube import clsim
    from os.path import expandvars, join, isfile
    from os import environ
    import tempfile, shutil
    import numpy
    
    propagator = I3CLSimStepToPhotonConverterPPC.instance
    propagator.SetGeometry(DetectorParams['Geometry'])
    propagator.SetWlenBias(DetectorParams['WavelengthGenerationBias'])

    tmpdir = tempfile.mkdtemp()

    # create config file
    with open(join(tmpdir, 'cfg.txt'), 'w') as f:
        write = lambda v: f.write('{}\n'.format(v))
        write(DetectorParams['DOMOversizeFactor'])
        write(DetectorParams['UnshadowedFraction']*DetectorParams['MediumProperties'].efficiency)

        # FIXME: extract scattering parameters from ice model
        icemodel = DetectorParams['MediumProperties']
        icemodel.GetScatteringCosAngleDistribution()
        write(0.45) # 0=HG; 1=SAM
        write(0.9) # g=<cos(theta)>
    
    shutil.copy(expandvars('$I3_BUILD/ppc/resources/ice/rnd.txt'), tmpdir)
    
    # configure unit angular acceptance, disabling angular downsampling in PPC
    with open(join(tmpdir, 'as.dat'), 'w') as f:
        for _ in range(2):
            f.write('{}\n'.format(1.0))

    # convert generation bias to cumulative distribution of photon wavelengths
    with open(join(tmpdir, 'wv.dat'), 'w') as f:
        # integrate Cherenkov spectrum weighted with given bias
        bias = DetectorParams['WavelengthGenerationBias']
        nPhase = DetectorParams['MediumProperties'].GetPhaseRefractiveIndex(0)
        def cherenkovDist(wlen, beta=1):
            # dN/dxdwlen
            return (2.*numpy.pi/(137.*(wlen*wlen)))*(1. - 1./ ( (beta*nPhase.GetValue(wlen))**2 ) )
        wavelength = numpy.array(list(map(bias.GetEntryWavelength, range(bias.GetNumEntries()))))
        value = numpy.array(list(map(cherenkovDist, wavelength)))*numpy.array(list(map(bias.GetEntryValue, range(bias.GetNumEntries()))))
        f.write('{:.8f} {:.0f}\n'.format(0, (wavelength[0] - bias.GetWavelengthStepping()/2)/I3Units.nanometer))
        for wvl, v in zip(wavelength + bias.GetWavelengthStepping()/2, numpy.cumsum(value)/sum(value)):
            f.write('{:.8f} {:.0f}\n'.format(v, wvl/I3Units.nanometer))
    
    # copy the ice model itself over blindly
    shutil.copy(join(DetectorParams['IceModelLocation'], 'icemodel.par'), tmpdir)
    shutil.copy(join(DetectorParams['IceModelLocation'], 'icemodel.dat'), tmpdir)
    
    for ext in 'par', 'dat':
        fname = join(DetectorParams['IceModelLocation'], 'tilt.'+ext)
        if isfile(fname):
            shutil.copy(fname, tmpdir)

    environ['PPCTABLESDIR'] = tmpdir
    if UseCPUs and not UseGPUs:
        environ['OCPU'] = '1'
    elif UseGPUs and not UseCPUs:
        environ['OGPU'] = '1'
    elif not UseGPUs and not UseCPUs:
        raise ValueError("Need to use at least one device type")
    
    try:
        propagator.Initialize()
    finally:
        shutil.rmtree(tmpdir)
    return propagator
