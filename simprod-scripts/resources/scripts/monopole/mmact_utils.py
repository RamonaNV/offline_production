from icecube import dataclasses, icetray
from icecube.icetray import I3Units
from icecube.dataclasses import I3MCTree, I3Particle

gen_params        = { "beta_spectrum":    [  0.750,
                                             0.995,  ],
                      "disk_distance":    1000.0    * I3Units.m,
                      "disk_radius":      1100.0    * I3Units.m,
                      "monopole_mass":       1.0e11 * I3Units.GeV,
                      "dist_to_cent_max": 2000.0    * I3Units.m,
                      "step_length_max":    10.0    * I3Units.m,   }

default_settings = { "n_events": 1000,
                     "icemodel": "$I3_BUILD/ice-models/resources/models/",
                     "gcd":      "/data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz", }

filename_template = "mmact_FLAVOR_IC86__beta_BETALOW_BETAHIGH__DATALEVEL_level__proc_PROCESSNUMBER.i3.gz"



def check_monopole_lengths_10m(frame):
	flavor = "Monopole"
	tree   = frame["I3MCTree"]
	maxlength = 10. * I3Units.m
	non_mm = [ int(str(mm.type)!=flavor) for mm in tree ]
	if any(non_mm):
		exit( "Found {} non-monopole particle(-s) in the MCTree! This should never happen!".format(sum(non_mm)) )
	too_long_mm = [ int(float(mm.length)>maxlength) for mm in tree ]
	if any(too_long_mm):
		exit( "Found {} too long monopole track segment(-s) in the MCTree! This is a bug that should have been fixed in the monopole-generator trunk!".format(sum(too_long_mm)) )

class mmact_print:
	def __init__(self,intime=-1,verb=1):
		self.script_verbose_level = verb
		self.label                = "MMACT"
		self.label_spaces         = "     "
		import time
		self.starttime            = intime if intime>0. else time.time()
		self.starttime_struct     = time.localtime(self.starttime)
		self.sts                  = self.starttime_struct
		self.endtime_struct       = -1
		self.ets                  = self.endtime_struct
	def set_label_spaces(self):
		self.label_spaces = " "*len(self.label)
	def start(self):
		self.set_label_spaces()
		self.sts = self.starttime_struct
		print(  " ----------------------------------------- " )
		print(  "| Starting script at  {yr:04}-{mo:02}-{dy:02} {h:02}:{m:02}:{s:02} |".format( yr=self.sts.tm_year, mo=self.sts.tm_mon, dy=self.sts.tm_mday, h=self.sts.tm_hour, m=self.sts.tm_min, s=self.sts.tm_sec ) )
		print(  " ----------------------------------------- " )
		print(  "" )
	def finish(self):
		self.set_label_spaces()
		import time
		self.endtime_struct = time.localtime(time.time())
		self.ets = self.endtime_struct
		print(  " -{l}- ----------  ".format( l=self.label_spaces.replace(" ","-") ) )
		print(  "" )
		print(  " ----------------------------------------- " )
		print(  "| Finishing script at {yr:04}-{mo:02}-{dy:02} {h:02}:{m:02}:{s:02} |".format( yr=self.ets.tm_year, mo=self.ets.tm_mon, dy=self.ets.tm_mday, h=self.ets.tm_hour, m=self.ets.tm_min, s=self.ets.tm_sec ) )
		print(  " ----------------------------------------- " )
	def vbprint(self,the_string,accepted_verbose_levels,info_level=0):
		if self.script_verbose_level in accepted_verbose_levels:
			self.set_label_spaces()
			self.label      = str(self.label)
			self.the_string = str(the_string)
			import time
			self.time_now    = int(time.time()-self.starttime)
			if info_level==0:
				print(  " -{l}- ----------  ".format(             l=self.label_spaces.replace(" ","-") ) )
				print(  "| {l} | {h:02}:{m:02}:{s:02} | ".format( l=self.label        , h=self.time_now/3600, m=(self.time_now%3600)/60, s=self.time_now%60 ) + the_string )
			if info_level==1:
				print(  "| {l} | {h:02}:{m:02}:{s:02} | ".format( l=self.label_spaces , h=self.time_now/3600, m=(self.time_now%3600)/60, s=self.time_now%60 ) + the_string )
			if info_level==2:
				print(  "| {l} | {h:02}:{m:02}:{s:02} | ".format( l=self.label_spaces , h=self.time_now/3600, m=(self.time_now%3600)/60, s=self.time_now%60 ) + the_string )


def add_MCPrimaryParticle(frame):
	if "I3MCTree" in frame and "MCPrimaryParticle" not in frame:
		accepted_primaries = ["monopole","nue","nuebar","numu","numubar","nutau","nutaubar"]
		try:
			primaries = [ p for p in frame["I3MCTree"].get_primaries() if str(p.type).lower() in accepted_primaries ]
		except AttributeError:
			primaries = [ p for p in [ frame["I3MCTree"][0] ] if str(p.type).lower() in accepted_primaries ]
		if len(primaries)!=1:
			exit("Event nr {} has {} primary particles (neutrinos and monopoles counted) in it's MCTree! Deal with this case!".format(frame["I3EventHeader"].event_id,len(primaries)))
		primary = primaries[0]
		frame["MCPrimaryParticle"] = primary


#
#def get_key_from_filename(template, filename, key):
#	filename = (filename.split("/")[-1]).replace(".","_")
#	template = (template               ).replace(".","_")
#	filename = filename.replace("i3_bz2","i3.bz2").replace("i3_gz","i3.gz")
#	temp_dict = { t: f for t, f in zip(template.split("_"),filename.split("_")) }
#	part = "YOUR_KEY_DOES_NOT_MATCH_ANY_KEY_IN_THE_TEMPLATE"
#	part = temp_dict[key.upper()]
#	return part
#
# def SOME_RANDOM_SEED
