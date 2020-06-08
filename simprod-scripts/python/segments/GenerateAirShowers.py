
from icecube.icetray import traysegment

@traysegment
def GenerateAirShowers(tray, name,
    Files=None,
    GCDFile=None,
    NEvents=1,
    OverSampling=1,
    LegacyOverSampling=False,
    SimulateIceTop=True,
    TankResponse="param",
    KeepMCHits=False,
    RandomService="I3RandomService",
    ):
	"""
	Read CORSIKA files, simulate IceTop response, and populate I3MCTree with penetrating components (neutrinos and muons)

	:param Files: list of CORSIKA files to read
	:param GCDFile: path to GCD file to read first
	:param NEvents: passed to I3CORSIKAReader (and ignored for CORSIKA files >= v74000, where it is part of the run header)
	:param OverSampling: Number of times each shower will be read (each with a different impact location)
	:param SimulateIceTop: simulate IceTop response
	:param TankResponse: how to simulate IceTop response: 'param' for parameterization, 'g4' for GEANT4
	:param KeepMCHits: keep IceTop MCHits
	:param RandomService: the name of a random service to be used by the tank response
	"""

	from icecube import icetray, corsika_reader
    	
	random = tray.context[RandomService]
	kwargs = dict()
	if SimulateIceTop:
		kwargs['ParticlesToWrite'] = [] # emit everything for IceTop simulation


	tray.Add(corsika_reader.ReadCorsika, "reader", 
		GCDFILE=GCDFile, FileNameList=Files, NEvents=NEvents, 
		OverSampling=OverSampling,LegacyOverSampling=LegacyOverSampling, **kwargs)


	if SimulateIceTop:
		from icecube import topsimulator
		# Simulate IceTop first
		if TankResponse.lower() == "g4":
			from icecube import g4_tankresponse
			tray.AddService("I3IceTopResponseFactory<I3G4TankResponse>", name+"topresponse",
					RandomServiceName = RandomService,
					ChargeScale=1.02
					)
		elif TankResponse.lower() == "param":
	        	tray.AddService("I3IceTopResponseFactory<I3ParamTankResponse>", name+"topresponse",
					RandomServiceName = RandomService,
					UseSnowParam=True
					)
		else:
			raise ValueError('Unknown tank response ' + str(TankResponse))
	
		tray.Add('topsim::MCTreeShimInjectorFactory', name+'topinjector')
	
	    	# The actual topsimulator module (depends on injector- and tank response services) 
		tray.AddModule("I3TopSimulator",name+"I3TopSimulator",
		    InjectorServiceName=name+"topinjector",
		    ResponseServiceName=name+"topresponse",
		    PrimaryName='', InIceMCTreeName='', IceTopMCTreeName='',
		)
	
		if not KeepMCHits:
			tray.Add("Delete", Keys=['IceTopMCHits'])


