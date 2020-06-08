
from icecube.icetray import traysegment, I3Units
from icecube import dataclasses

@traysegment
def GenerateBundles(tray, name, Generator=None,
    RunNumber=1, NEvents=100,
    GCDFile='/data/sim/sim-new/downloads/GCD_31_08_11/GeoCalibDetectorStatus_IC79.55380_L2a.i3.gz',
    FromTime=dataclasses.I3Time(55380),
    ToTime=dataclasses.I3Time(55380)):
	"""
	Generate muon bundles from a parametrization.
	
	:param Generator: an instance of I3MuonGun.Generator
	"""
	
	from icecube import icetray, dataclasses
	from icecube import sim_services, MuonGun
	
	RandomService = tray.context['I3RandomService']
	
	tray.AddModule("I3InfiniteSource",name+"_streams",
	    Prefix=GCDFile, Stream=icetray.I3Frame.DAQ)
		
	# modify the header if necessary by generating a random time
	# between StartTime and EndTime
	def GenRandomTime(frame, StartTime, EndTime, RandomService):
		header = dataclasses.I3EventHeader()

		mjdStart = StartTime.mod_julian_day_double
		mjdEnd	 = EndTime.mod_julian_day_double
		
		mjd = mjdStart + (RandomService.uniform(0.,1.)*(mjdEnd-mjdStart))
		eventTime = dataclasses.I3Time()
		eventTime.set_mod_julian_time_double(mjd)
		header.start_time = eventTime
		frame["I3EventHeader"] = header
		
	if FromTime != ToTime:
		if RandomService is None:
			raise ValueError("You must specify a random service!")
		tray.AddModule(GenRandomTime, name+"_GenRandomTime",
			StartTime=FromTime,
			EndTime=ToTime,
			RandomService=RandomService,
			Streams=[icetray.I3Frame.DAQ])
	
	tray.AddModule('I3MuonGun::GeneratorModule',name,
	    Generator=NEvents*Generator)
	# tray.AddModule('Dump', name+'dump')
	
