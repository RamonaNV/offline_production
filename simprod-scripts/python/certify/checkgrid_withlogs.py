#!/usr/bin/env python
################################ Tray 0 ##############################
#   on 2010-01-11 
#   
######################################################################
import commands,os,sys,string
from os.path import expandvars
import glob
import logging

from I3Tray import *
from iceprod.modules import ipmodule

# Instantiate parameters expression parser 

def boolcast(s): return s not in ["False","0","","[]"]
double = float

class CertifyGrid(ipmodule.IPBaseClass):
   """
   SysChk for production grids
   """
   def Execute(self,stats):
	if not ipmodule.IPBaseClass.Execute(self,stats):
           return 0
	
        import icecube.icetray
	import cPickle
	import random
	import md5 
	return 0

class CopyLogs(ipmodule.IPBaseClass):
   """
   Copy log/outfiles to storage units with GridFTP
   """
   def Execute(self,stats):
 	if not ipmodule.IPBaseClass.Execute(self,stats): return 0
	
	import icecube.icetray
        from iceprod.modules import gsiftp
        from iceprod.core.dataclasses import I3PostTray

        post    = I3PostTray()
        job     = int(self.parser.parse("$args(procnum)"))
        key     = self.parser.parse("$args(key)")
        dataset = int(self.parser.parse("$args(dataset)"))
	TopDir = expandvars("$TMPTOP")	
	
	station = self.parser.parse("$system(gridname)")
	tx_protocol = "gsiftp://data1.icecube.wisc.edu"
	if str(station) == "npx3": tx_protocol = "file:"
				
	os.system(""" tar -cvf tarred_file.tar %s/*log* %s/*out* %s/*err*"""%(TopDir,TopDir,TopDir))	
		
	post.AddModule("gsiftp.URLCopy","upload_tar")(
	("source","file:%s" % os.path.join(os.getcwd(),"tarred_file.tar")),
    	("destination","%s/data/user/ice3simusr/SysChkLogs/%s/%s/%s_%s_tarred_file.tar"%(tx_protocol,str(station),str(dataset),str(job),str(key))),
	)
		
	post.Execute()

	logging.debug("*******")
	logging.debug(job)
	logging.debug(key)
	logging.debug("*********")
	raise ("!")

	return 0
