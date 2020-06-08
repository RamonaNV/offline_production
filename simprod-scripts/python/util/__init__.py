from icecube import icetray


def DrivingTime( frame ):
    time = frame["I3EventHeader"].start_time
    if "DrivingTime" not in frame : 
       frame.Put("DrivingTime", time )

def AddEmptyWeightMap(frame, mapname = "CorsikaWeightMap" ):
    from icecube import dataclasses
    weights = dataclasses.I3MapStringDouble()
    weights["Weight"] = 1.0
    if mapname not in frame:
       frame.Put(mapname,weights)
       


def CombineHits(frame,
        InputHits = ["I3MCPESeriesMapTracks", "I3MCPESeriesMapCascades"],
        OutputHits = "I3MCPESeriesMap"):
    from icecube import dataclasses, polyplopia, simclasses

    newhits = simclasses.I3MCPESeriesMap()
    for hitmap in InputHits:
        hits = frame[hitmap]
        polyplopia.MergeHits(newhits,hits,0.0)
    frame[OutputHits] = newhits
    return True
         


def BasicHitFilter(frame,HitSeriesMapNames="I3MCPESeriesMap"):
    if frame.Has(HitSeriesMapNames):
       hits = frame.Get(HitSeriesMapNames)
       return sum(map(len,hits.values()))
    return False



def BasicCounter(frame,name="basic_counter",Stats={}):
      if name not in Stats: Stats[name]=0
      Stats[name] += 1

class DAQCounter(icetray.I3Module):
    count   = 0
    nevents = 0

    def __init__(self,ctx):
        icetray.I3Module.__init__(self,ctx)
        self.AddParameter("NEvents","name of event counter",self.nevents)
        self.AddOutBox("OutBox");

    def Configure(self):
        self.nevents = self.GetParameter("NEvents")

    def DAQ(self,frame):
        self.count += 1
        if self.count > self.nevents: self.RequestSuspension()
        self.PushFrame(frame)


class PrintContext(icetray.I3Module):

    def __init__(self,ctx):
        icetray.I3Module.__init__(self,ctx)
        self.ctx = ctx
        self.AddOutBox("OutBox");

    def Configure(self):
        icetray.logging.log_info("Keys: {}".format(self.ctx.keys()),"Context")



def choose_max_efficiency(efficiencies):
    """
    Given a single or list of DOM efficiencies choose the highest
    """
    if type(efficiencies) == list or type(efficiencies) == tuple:
       return max(map(float,efficiencies)) 
    else:
       return float(efficiencies)
 


def I3Summary2JSON(summary):
    from icecube.dataclasses import I3MapStringDouble
    import json
    return json.dumps(dict(zip(summary.keys(),summary.values())))


def JSON2I3Summary(j):
    from icecube.dataclasses import I3MapStringDouble
    import json
    d = json.loads(j)
    summary =  I3MapStringDouble()
    for key,value in d.items():
        summary[str(key)] = value
    return summary



def WriteI3Summary(summary, filename):
    from icecube.dataclasses import I3MapStringDouble
    import json
    fs = open(filename,'w')
    return json.dump(dict(zip(summary.keys(),summary.values())), fs)

def ReadI3Summary(filename):
    from icecube.dataclasses import I3MapStringDouble
    import json
    fs = open(filename,'r')
    return JSON2I3Summary(fs.read())

def I3Summary2Dict(s):
    from icecube.dataclasses import I3MapStringDouble
    return dict(zip(summary.keys(),summary.values()))


def SetGPUEnvironmentVariables(GPU):
    os.putenv("CUDA_VISIBLE_DEVICES", str(GPU))
    os.putenv("COMPUTE", ":0." + str(GPU))
    os.putenv("GPU_DEVICE_ORDINAL", str(GPU))

