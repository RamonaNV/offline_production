"""
 * class: ModifyTime
 * (c) 2014 IceCube Collaboration
 * author: Juan Carlos Diaz Velez (juancarlos@icecube.wisc.edu)
 * IceProd classes for modifiying event frames
 * 
"""
######################################################################
import os
import sys,getopt
from os.path import expandvars

from icecube import icetray,dataclasses
from ..util import ReadI3Summary, WriteI3Summary
from .. import ipmodule

class ModifyEventTime(icetray.I3Module):
   """
   Modify the time in the I3EventHeader
   """
   def __init__(self,ctx) : 
       icetray.I3Module.__init__(self,ctx)

       self.AddParameter("MJD","Modified Julian Date",55697);
       self.AddParameter("MJDSeconds","Number of seconds after the start of the MJD.",0);
       self.AddParameter("MJDNanoSeconds","Number of nanoseconds after the start of the MJD.",0);
       self.AddParameter("Randomize","Randomize mjd within year starting with MJD",False);
       self.AddParameter("RNG","(pseudo-) random number generator",None);
       self.AddOutBox("OutBox");

   def Configure(self):
       mjd     = self.GetParameter("MJD")
       sec     = self.GetParameter("MJDSeconds")
       nanosec = self.GetParameter("MJDNanoSeconds")
       self.rand = self.GetParameter("Randomize") 
       self.rng  = self.GetParameter("RNG") 
       self.time  = dataclasses.I3Time()
       self.time.set_mod_julian_time(mjd,sec,nanosec)

   def DAQ(self,frame):
       if self.rand and self.rng is not None:
          mjdStart = self.GetParameter("MJD")
          mjdEnd = self.GetParameter("MJD") + 365
          mjd = self.rng.uniform(0.,1.)*(mjdEnd-mjdStart)
          sec = self.rng.uniform(0.,1.)*24*3600
          self.time.set_mod_julian_time(mjd,sec,0)

       if "I3EventHeader" in frame:
          eventheader = frame["I3EventHeader"]
          eventheader.start_time = self.time
          del frame["I3EventHeader"]
          frame["I3EventHeader"] = eventheader
       else:
          icetray.logging.log_warn("No I3EventHeader to modify.")

       # some GCDs come with driving times and some don't
       if "DrivingTime" in frame :
          del frame["DrivingTime"]
       frame.Put("DrivingTime", self.time )
       
       self.PushFrame(frame,"OutBox")

 

class DateShifter(ipmodule.ParsingModule):

   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('inputfile','Output filename','$steering(prev_file)')
        self.AddParameter('outputfile','Output filename','$steering(current_file)')
        self.AddParameter('summaryfile','JSON Summary filename','summary.json')
        self.AddParameter("MJD","Modified Julian Date",55697)
        self.AddParameter("MJDSeconds","Time added in seconds after mjd",0)
        self.AddParameter("MJDNanoSeconds","Time added in seconds after mjd",0)

   def Execute(self,stats):
        if not ipmodule.IPBaseClass.Execute(self,stats): return 0

        import icecube.icetray
        from I3Tray import I3Tray
        from icecube import dataio,dataclasses,phys_services

        inputfile = self.GetParameter('inputfile')
        outputfile = self.GetParameter('outputfile')


        # Instantiate a tray 
        tray = I3Tray()

        # Configure IceTray services
        summary = dataclasses.I3MapStringDouble()
        summaryfile     = self.GetParameter('summaryfile')
        if os.path.exists(summaryfile): 
           summary = ReadI3Summary(summaryfile)
        tray.context['I3SummaryService'] = summary
         
        # Configure IceTray modules 
        tray.AddModule("I3Reader","reader",filenamelist=[inputfile])
        tray.AddModule(ModifyEventTime,"mod_mjd", 
             MJD=self.GetParameter('mjd'),
             MJDSeconds=self.GetParameter('mjdseconds'),
             MJDNanoSeconds=self.GetParameter('mjdnanoseconds'))
        tray.AddModule("I3Writer","writer", filename=outputfile)
        

        # Execute the Tray
        tray.Execute()
        

        tray.PrintUsage()

        summary = tray.context['I3SummaryService']

        # Save stats
        for k in tray.Usage():
            stats[str(k.key())+":usr"] = k.data().usertime
            summary[str(k.key())+":usr"] = k.data().usertime

            stats[str(k.key())+":sys"] = k.data().systime
            summary[str(k.key())+":sys"] = k.data().systime

            stats[str(k.key())+":ncall"] = k.data().ncall
            summary[str(k.key())+":ncall"] = k.data().ncall

        WriteI3Summary(summary, summaryfile)

        # Free memory
        del tray
        return 0

