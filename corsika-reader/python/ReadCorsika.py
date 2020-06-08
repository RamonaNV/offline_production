
from icecube.icetray import traysegment

@traysegment
def ReadCorsika(tray, name,
    filenamelist=[],
    GCDFile='',
    NEvents=1,
    OverSampling=1,
    LegacyOverSampling=False,
    CylinderHeight=1200,
    CylinderRadius=600,
    RandomService="I3RandomService",
    TrimShower=True
    ):
	"""
	Read CORSIKA files, simulate IceTop response, and populate I3MCTree with penetrating components (neutrinos and muons)

	:param Files: list of CORSIKA files to read
	:param GCDFile: path to GCD file to read first
	:param NEvents: passed to I3CORSIKAReader (and ignored for CORSIKA files >= v74000, where it is part of the run header)
	:param OverSampling: Number of times each shower will be read (each with a different impact location)
	:param CylinderHeight: height of IceCube-centered target cylinder
	:param CylinderRadius: radius of target cylinder
	:param RandomService: the name of a random service to be used by the tank response
	:param TrimShower: remove surface particles from tree
	"""

	from icecube import icetray, corsika_reader, sim_services
    	
	random = tray.context[RandomService]
	kwargs = dict()

	legacy_ov_samp = 1
	ov_samp = OverSampling
	if LegacyOverSampling: 
		legacy_ov_samp = OverSampling
		ov_samp = 1

	tray.Add('I3CORSIKAReader', 'reader', filenamelist=filenamelist, NEvents=NEvents,
	    Prefix=GCDFile, LegacyOverSampling=legacy_ov_samp, **kwargs)

	if TrimShower:
		# Then drop all the surface crap
		tray.Add('I3InIceCORSIKATrimmer')
	# Drop showers where no particles reach the observation level
	tray.Add(lambda frame: len(frame['I3MCTree']) > 1, Streams=[icetray.I3Frame.DAQ])

	tray.Add('Rename','rename_tree',keys=['I3MCTree','I3MCTree_preSampling'])

	if OverSampling > 0:
		tray.AddModule("CORSIKAResampler","resample",
	        	OverSampling=ov_samp,
	        	CylinderHeight=CylinderHeight, CylinderRadius=CylinderRadius)



