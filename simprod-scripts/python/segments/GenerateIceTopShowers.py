from icecube.icetray import traysegment, I3Units

@traysegment
def GenerateIceTopShowers(tray, name,
			  Files=None,
			  GCDFile=None,
			  NSamples=1,
			  TankResponse="g4",
			  KeepMCHits=False,
			  x=0, y=0, r=0,
			  ResponseRandomService="I3RandomService",
			  InjectorRandomService="I3RandomService",
			  UnthinRadius=0,
			  TankSamplingRadius=25*I3Units.m,
			  RunID=0,
			  RaiseObservationLevel=0
			  ):
	"""
	Read CORSIKA files, simulate IceTop response, and populate I3MCTree with penetrating components (neutrinos and muons) using topsimulator's injector.

	:param Files: list of CORSIKA files to read
	:param GCDFile: path to GCD file to read first
	:param NSampless: Number of samples of a single CORSIKA shower
	:param TankResponse: Method used to simulate IceTop response: 'param' for parameterization, 'g4' for GEANT4
	:param KeepMCHits: keep IceTop MCHits (needed only for the PMTSimulator-based simulation)
	:param x: x-coordinate of the point around which the CORSIKA origin is randomly placed
	:param y: y-coordinate of the point around which the CORSIKA origin is randomly placed
	:param r: radius of the circle in which the CORSIKA origin is randomly placed
	:param ResponseRandomService: name of the I3RandomService used for tank response
	:param InjectorRandomService: name of the I3RandomService used for particle injection
	:param UnthinRadius: radius of unthinning sampling area
	:param TankSamplingRadius: the radius around tanks where all particles are considered for propagation (>=UnthinRadius)
	:param RunID: the ID to store in the header (normally that of the corsika run)
	"""

	from icecube import icetray, dataio
	from icecube import topsimulator
	import os

	# we compres CORSIKA files now. We need to change our reader to handle this. Meanwhile...
	# this assumes we are copying the file before, for this reason I added a line to check the relative path
	missing = []
	for f in Files:
		if not os.path.dirname(os.path.relpath(f)):
			if os.path.exists('%s.bz2'%f):
				rc = os.system("bzip2 -d '" + f + ".bz2' >&2")
			if os.path.exists('%s.gz'%f):
				rc = os.system("gzip -d '" + f + ".gz' >&2")
		if not os.path.exists(f):
			missing.append(f)
	if missing:
		raise Exception("Missing input files: %s"%str(missing))

	tray.AddModule("I3InfiniteSource", "source",
		       prefix = GCDFile,
		       stream = icetray.I3Frame.DAQ )

	from icecube.sim_services.sim_utils.gcd_utils import get_time
	time = get_time(dataio.I3File(GCDFile))

	# topsimulator overrides the given I3EventHeader with this time... (?)
	def DrivingTime( frame, t=time ):
		if "DrivingTime" in frame :
			del frame["DrivingTime"]
		frame.Put("DrivingTime", t )

	tray.AddModule(DrivingTime, "dt",  
		       Streams = [icetray.I3Frame.DAQ] ) 

	tray.AddModule("I3MCEventHeaderGenerator","time-gen",
		Year = time.utc_year,
		DAQTime = time.utc_daq_time,
		RunNumber = RunID,
		IncrementEventID = True
		)

	if TankResponse.lower() == "g4":
		from icecube import g4_tankresponse
		tray.AddService("I3IceTopResponseFactory<I3G4TankResponse>", name+"_topresponse",
				RandomServiceName = ResponseRandomService,
				ChargeScale=1.02, # From Arne. Javier doesn't believe this is exactly right...
				TankSamplingRadius=TankSamplingRadius,
				)
	elif TankResponse.lower() == 'param':
		tray.AddService("I3IceTopResponseFactory<I3ParamTankResponse>", name+"_topresponse",
				RandomServiceName = ResponseRandomService,
				UseSnowParam = True,
				TankSamplingRadius=TankSamplingRadius,
				)
	else:
		raise Exception('Unknown tank response ' + str(TankResponse))
	
	tray.AddService("I3InjectorFactory<I3CorsikaInjector>", name+"_topinjector",
			FileNameList = Files,
			RandomServiceName = InjectorRandomService,
			NumSamples = NSamples,
			RelocationX = x,
			RelocationY = y,
			RelocationR = r,
			UnThinRadius = UnthinRadius,
			IgnoreParticleTypes = [75, 76, 85, 86, 95, 96],
			RaiseObservationLevel = RaiseObservationLevel
			)

	tray.AddModule("I3TopSimulator",name+"I3TopSimulator",
		       InIceMCTreeName="I3MCTree",
		       IceTopMCTreeName="IceTopMCTree",
		       ResponseServiceName=name+"_topresponse",
		       InjectorServiceName=name+"_topinjector",
		       PrimaryName='MCPrimary',
		       IceTopHitSeriesName='IceTopMCHits',
		       MuonEnergyCutOff=0., # this is done by I3InIceCORSIKATrimmer, so it can be anything.
		       CompressPEs=1,
		       UseInjectorComponents=True
		       )
	
	tray.AddModule("I3EventCounter","top_simulator_counter",
		       physicscountername="Generated Events",
		       CounterStep=1,
		       )

	if not KeepMCHits:
		tray.Add("Delete", Keys=['IceTopMCHits'])

	tray.Add('I3InIceCORSIKATrimmer', FilterMode=False)

