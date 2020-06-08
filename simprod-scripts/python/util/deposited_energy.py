
import os
from icecube import icetray, dataclasses, simclasses

class I3DepositedEnergy(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        #self.context = context
        self.AddParameter("MMCListName", "Name of I3MMCList object", "I3MMCList")
	    
    def Configure(self):
        self.mmclistname = self.GetParameter("MMCListName")
        self.emax = 0.

    def DAQ(self, frame):
        
        mmclist = frame[self.mmclistname]
        etot = 0.
        for mmc in mmclist:
        	etot += mmc.Elost
        print etot
        self.emax = max(etot,self.emax)
        self.PushFrame(frame)

    def Finish(self):
        summary = self.context['I3SummaryService'] 
        summary['MaxEnergyDeposited'] = self.emax
        print("Maximum energy deposited: %f" % self.emax )

