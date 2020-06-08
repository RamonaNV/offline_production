#!/usr/bin/env python

from os.path import expandvars
from icecube import icetray, dataclasses, phys_services
from I3Tray import I3Units

import math

@icetray.traysegment
def GenerateNeutrinos(tray, name,
    RandomService=None,
    RunID=None,
    NumEvents=100,                
    SimMode='Full',                
    VTXGenMode='NuGen',         # currently only NuGen is supported
    InjectionMode='Surface',
    CylinderParams=[0,0,0,0,0], # CIRCLE[radius, active_height_before, active_height_after] SURFACE[radius, length, center_x, center_y, center_z]
    AutoExtendMuonVolume=True,   
    Flavor="", # if it is set, I3NuGInjector is used.
    NuTypes=["NuMu","NuMuBar"], # if it is set, I3NuGDiffuseSource is used.              
    PrimaryTypeRatio=[1,1],               
    GammaIndex=2.0,             
    FromEnergy=1.*I3Units.TeV,
    ToEnergy=10.*I3Units.PeV,
    ZenithRange=[0, 180*I3Units.degree],
    AzimuthRange=[0, 360*I3Units.degree],
    UseDifferentialXsection=True,           # set true if you use differential cross section
    CrossSections='csms_differential_v1.0', # 
    CrossSectionsPath=None, 
    ZenithSamplingMode='ANGEMU', 
    ParamsMap = dict() 
    ):

         
    #-----------------------------
    # parameter check
    #-----------------------------

    if Flavor == "" and len(NuTypes) == 0 :
        raise Exception('Need to set Flavor or NuTypes.')

    #-----------------------------
    # icetray start
    #-----------------------------

    from I3Tray import I3Units

    from icecube import icetray, dataclasses, dataio, phys_services, sim_services
    from icecube import neutrino_generator, MuonGun
    
    earthModelService = name+"_earthmodel"
    steering = name+"_steering"
    injector = name+"_injector"
    interactions = name+"_interactions"

    #-----------------------------
    # make ParamMap's key lowercase and
    # check multiple definitions.
    #-----------------------------

    params = dict()
    for key in ParamsMap:
        keylower = key.lower()
        if keylower in params:
            raise Exception('param %s is set twice!'%(keylower))
        params[keylower] = ParamsMap[key]

    #-----------------------------
    # configure earthmodel service 
    #-----------------------------

    earthModelServiceArgs = dict()
    if "earthmodels" in params:
        if not isinstance(params["earthmodels"],(list,tuple)):
            raise Exception('EarthModels must be a list of strings')
        earthModelServiceArgs['EarthModels'] = params["earthmodels"]

    if "materialmodels" in params :
        if not isinstance(params["materialmodels"],(list,tuple)):
            raise Exception('MaterialModels must be a list of strings')
        earthModelServiceArgs['MaterialModels'] = params["materialmodels"] 

    if "earthparamspath" in params:
        params["earthparamspath"] = str(params["earthparamspath"])
        earthModelServiceArgs['PathToDataFileDir'] = params["earthparamspath"] 

    if "icecaptype" in params:
        params["icecaptype"] = str(params["icecaptype"])
        earthModelServiceArgs['IceCapType'] = params["icecaptype"] 

    if "icecapsimpleangle" in params:
        if not isinstance(params["icecapsimpleangle"], float):
            raise Exception('IceCapType must be a float')
        earthModelServiceArgs['IceCapSimpleAngle'] = params["icecapsimpleangle"]

    if "detectordepth" in params :
        if not isinstance(params["detectordepth"], float):
            raise Exception('DetectorDepth must be a float')
        earthModelServiceArgs['DetectorDepth'] = params["detectordepth"]

    tray.AddService("I3EarthModelServiceFactory", earthModelService,
            **earthModelServiceArgs
        )

    #-----------------------------
    # configure steering factory
    #-----------------------------

    # this is default cylinder. If you want to change cylinder size
    # and cylinder center, set positive value to CylinderParams.
    surface = phys_services.Cylinder(1900*I3Units.m, 950*I3Units.m)
    #surface = MuonGun.Cylinder(1900*I3Units.m, 950*I3Units.m)

    steeringServiceArgs = dict()
    steeringServiceArgs['MCTreeName'] = 'I3MCTree_preMuonProp'
    steeringServiceArgs['NEvents'] = NumEvents
    steeringServiceArgs['SimMode'] = SimMode 
    steeringServiceArgs['VTXGenMode'] = VTXGenMode 
    steeringServiceArgs['InjectionMode'] = InjectionMode
    steeringServiceArgs['DetectionSurface'] = surface
    steeringServiceArgs['CylinderParams'] = CylinderParams 
    steeringServiceArgs['DoMuonRangeExtension'] = AutoExtendMuonVolume

    if "globalxsecscalefactor" in params :
        if not isinstance(params["globalxsecscalefactor"],(list,tuple)):
            raise Exception('GlobalXsecScaleFactor must be a list of float')
        steeringServiceArgs['GlobalXsecScaleFactor'] = params["globalxsecscalefactor"] 

    if "usesimplescatterform" in params :
        if not isinstance(params['usesimplescatterform'], int):
            raise Exception('UseSimpleScatterForm must be an int ')
        if params['usesimplescatterform'] > 0 :
            steeringServiceArgs['UseSimpleScatterForm'] = True
        else :
            steeringServiceArgs['UseSimpleScatterForm'] = False
 
    tray.AddService("I3NuGSteeringFactory", steering,
        EarthModelName=earthModelService,
        **steeringServiceArgs
        )


    #-----------------------------
    # configure injector
    #-----------------------------
    
    EnergyLogRange=[math.log10(FromEnergy/I3Units.GeV),math.log10(ToEnergy/I3Units.GeV)]    
    
    injectorServiceArgs = dict()
    injectorServiceArgs['GammaIndex'] = GammaIndex
    injectorServiceArgs['ZenithMin'] = ZenithRange[0]
    injectorServiceArgs['ZenithMax'] = ZenithRange[1]
    injectorServiceArgs['AzimuthMin'] = AzimuthRange[0]
    injectorServiceArgs['AzimuthMax'] = AzimuthRange[1]
    injectorServiceArgs['EnergyMinLog'] = EnergyLogRange[0]
    injectorServiceArgs['EnergyMaxLog'] = EnergyLogRange[1]
    injectorServiceArgs['AngleSamplingMode'] = ZenithSamplingMode
    injectorServiceArgs['RandomService'] = RandomService

    if "zenithweightparam" in params :
        if not isinstance(params["zenithweightparam"], float):
            raise Exception('ZenithWeightParam must be a float')
        injectorServiceArgs['ZenithWeightParam'] = params["zenithweightparam"] 

    if Flavor != "" :
        injectorServiceArgs['NuFlavor'] = Flavor
        tray.AddService("I3NuGInjectorFactory", injector,
            SteeringName=steering,
            **injectorServiceArgs
            )
    else :
        injectorServiceArgs['NuTypes'] = NuTypes
        injectorServiceArgs['PrimaryTypeRatio'] = PrimaryTypeRatio 
        tray.AddModule("I3NuGDiffuseSource", injector,
            SteeringName=steering,
            **injectorServiceArgs
            )

    #-----------------------------
    # configure interaction service
    #-----------------------------

    interactionServiceArgs = dict()
    if CrossSectionsPath is not None :
        interactionServiceArgs['TablesDir'] = CrossSectionsPath
    interactionServiceArgs['CrossSectionModel'] = CrossSections

    if UseDifferentialXsection == True :
        tray.AddService("I3NuGInteractionInfoDifferentialFactory", interactions,
            SteeringName=steering,
            **interactionServiceArgs
            )
    else :
        tray.AddService("I3NuGInteractionInfoFactory", interactions,
            SteeringName=steering,
            **interactionServiceArgs
            )

    #-----------------------------
    # configure neutrino generator
    #-----------------------------

    nugenArgs = dict()

    if "primarynuname" in params :
        params["primarynuname"] = str(params["primarynuname"])
        nugenArgs['PrimaryNuName'] = params["primarynuname"]

    if "interactionweight" in params :
        if not isinstance(params["interactionweight"],(list,tuple)):
            raise Exception('InteractionWeight must be a list of float')
        nugenArgs['InteractionCCFactor'] = params["interactionweight"][0] 
        nugenArgs['InteractionNCFactor'] = params["interactionweight"][1] 
        nugenArgs['InteractionGRFactor'] = params["interactionweight"][2] 
        nugenArgs['RandomService'] = RandomService

    if "propagationweightmode" in params :
        params["propagationweightmode"] = str(params["propagationweightmode"])
        propmode = neutrino_generator.to_propagation_mode(params["propagationweightmode"])
        nugenArgs['PropagationWeightMode'] = propmode

    if RunID is not None:
        nugenArgs['RunID']=RunID
        
    tray.AddModule("I3NeutrinoGenerator",name+"_neutrino",
        SteeringName=steering,
        InjectorName=injector,
        InteractionInfoName=interactions,
        **nugenArgs
        )


@icetray.traysegment
def SelectNeutrino(tray, name,
    Propagators=None,
    AutoExtendMuonVolume=False,
    EnergyBiasPower=1,
    FlavorBias=[30,1,1],
    CylinderRadius=600*I3Units.m,
    CylinderHeight=1200*I3Units.m,
    CrossSections='csms',
    ):
    r"""
    Select a neutrino to interact, and add neutrino propagators to the context.

    :param AutoExtendMuonVolume: allow :math:`\nu_{\mu}` to interact far before they reach the detector
    :param EnergyBiasPower: select a neutrino from the bundle with probability proportional to E^power
    :param FlavorBias: scale selection probability for :math:`\nu_{e}/\nu_{\mu}/\nu_{\tau}`
                       by these factors. The default value is appropriate for equal sampling
                       of conventional atmospheric :math:`\nu_{e}/\nu_{\mu}`.
    :param CylinderRadius: radius of upright-cylinder target volume
    :param CylinderHeight: full height of simulation volume
    :param CrossSections: cross-section tables to use ('cteq5', 'css', or 'csms')
    """
    from operator import add
    from icecube import neutrino_generator, sim_services, MuonGun

    # Set up NeutrinoGenerator internals
    random = tray.context['I3RandomService']
    #surface = MuonGun.Cylinder(CylinderHeight, CylinderRadius)
    surface = phys_services.Cylinder(CylinderHeight, CylinderRadius)
    config = neutrino_generator.Steering()
    config.detection_surface = surface
    config.do_muon_range_extension = AutoExtendMuonVolume
    interactions = neutrino_generator.I3NuGInteractionInfo(random, config, CrossSections)
    interactions.initialize()
    # I3NuGSourceSelector needs this in its context
    tray.context['NuGSteer'] = config
    
    # Remove all but one neutrino
    tray.Add('I3NuGSourceSelector', EnergyBiasPowerIndex=EnergyBiasPower,
        ParticleBiases=reduce(add, [[b]*2 for b in FlavorBias]),
        KeepDarkNeutrinos=False)

    # Store propagators in the context
    if not 'I3ParticleTypePropagatorServiceMap' in tray.context:
        tray.context['I3ParticleTypePropagatorServiceMap'] = sim_services.I3ParticleTypePropagatorServiceMap()
    Propagators = tray.context['I3ParticleTypePropagatorServiceMap']

    # Use NeutrinoPropagator for neutrinos
    prop = neutrino_generator.I3NeutrinoPropagator(random, config, interactions)
    # ensure that all nu_e and nu_mu reach the detector
    prop.prop_mode = neutrino_generator.PropagationMode.ncgrweighted
    tau_prop = neutrino_generator.I3NeutrinoPropagator(random, config, interactions)
    # un-weighted propagation for nu_tau to allow for tau regeneration
    tau_prop.prop_mode = neutrino_generator.PropagationMode.nopropweight
    for flavor in 'E', 'Mu', 'Tau':
        for ptype in '', 'Bar':
            Propagators[getattr(dataclasses.I3Particle.ParticleType, 'Nu'+flavor+ptype)] = tau_prop if flavor == 'Tau' else prop

@icetray.traysegment
def GenerateAtmosphericNeutrinos(tray, name,
    Files,
    GCDFile="",
    AutoExtendMuonVolume=False,
    EnergyBiasPower=1,
    FlavorBias=[30,1,1],
    CylinderRadius=600*I3Units.m,
    CylinderHeight=1200*I3Units.m,
    CrossSections='csms',
    NEvents=-1,
    MakeNeutrino=True
    ):
    r"""
    Read CORSIKA showers containing neutrinos, and force exactly one neutrino to interact.

    NB: this segment is deprecated. Use GenerateAirShowers, SelectNeutrino, and PropagateMuons instead.
    
    :param Files: a list of CORSIKA files to read
    :param GCDFile: I3 file with GCD information to read in before the CORSIKA files
    :param AutoExtendMuonVolume: allow :math:`\nu_{\mu}` to interact far before they reach the detector
    :param EnergyBiasPower: select a neutrino from the bundle with probability proportional to E^power
    :param FlavorBias: scale selection probability for :math:`\nu_{e}/\nu_{\mu}/\nu_{\tau}`
                       by these factors. The default value is appropriate for equal sampling
                       of conventional atmospheric :math:`\nu_{e}/\nu_{\mu}`.
    :param CylinderRadius: radius of upright-cylinder target volume
    :param CylinderHeight: full height of simulation volume
    :param CrossSections: cross-section tables to use ('cteq5', 'css', or 'csms')
    """
    import warnings
    warnings.warn('GenerateAtmosphericNeutrinos is deprecated. Use GenerateAirShowers, SelectNeutrino, and PropagateMuons instead')
    from operator import add
    from icecube import neutrino_generator, sim_services, MuonGun
    from icecube.sim_services.propagation import get_propagators
    icetray.load('corsika-reader', False)

    random = tray.context['I3RandomService']
    surface = MuonGun.Cylinder(CylinderHeight, CylinderRadius)
    tray.Add('I3CORSIKAReader', 'reader', filenamelist=Files, NEvents=NEvents,
        CylinderHeight=surface.length, CylinderRadius=surface.radius, Prefix=GCDFile)

    # Drop showers where no particles reach the observation level
    tray.Add(lambda frame: len(frame['I3MCTree']) > 1, Streams=[icetray.I3Frame.DAQ])

    # Remove low-energy muons that can't reach the detector
    tray.Add('I3InIceCORSIKATrimmer')

    tray.Add(SelectNeutrino, AutoExtendMuonVolume=AutoExtendMuonVolume, EnergyBiasPower=EnergyBiasPower,
        FlavorBias=FlavorBias, CylinderRadius=CylinderRadius, CylinderHeight=CylinderHeight,
        CrossSections=CrossSections)
    
    base_propagators = get_propagators()
    propagators = tray.context['I3ParticleTypePropagatorServiceMap']
    for k in base_propagators.keys():
        propagators[k] = base_propagators[k]
    
    tray.Add('Rename', Keys=['I3MCTree', 'I3MCTree_preMuonProp'])
    tray.Add('I3PropagatorModule', PropagatorServices=propagators,
             InputMCTreeName="I3MCTree_preMuonProp", OutputMCTreeName="I3MCTree",
             RandomService='I3RandomService')
