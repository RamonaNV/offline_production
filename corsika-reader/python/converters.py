#print(__name__,__path__)
from icecube import tableio
#from . import I3CorsikaWeight





#from icecube.load_pybindings import load_pybindings

from icecube.corsika_reader import I3CorsikaWeight
from icecube.dataclasses import I3Particle

class I3CorsikaWeightConverter(tableio.I3Converter):
    booked = I3CorsikaWeight
    def CreateDescription(self,lfparams):
        desc = tableio.I3TableRowDescription()
        desc.add_field("primary_zenith",tableio.types.Float64,"","")
        desc.add_field("primary_azimuth",tableio.types.Float64,"","")
        desc.add_field("primary_energy",tableio.types.Float64,"","")
        desc.add_field('primary_type', tableio.I3Datatype(I3Particle.ParticleType), '', '')
        desc.add_field("bias_factor",tableio.types.Float64,"","")
        desc.add_field("bias_target",tableio.types.Float64,"","")
        desc.add_field("weight",tableio.types.Float64,"","")
        desc.add_field("max_x",tableio.types.Float64,"","")                
        return desc
    def FillRows(self,lfparams,rows):

        rows["primary_zenith"]  = lfparams.primary.dir.zenith
        rows["primary_azimuth"] = lfparams.primary.dir.azimuth
        rows["primary_energy"]  = lfparams.primary.energy
        rows["primary_type"]    = lfparams.primary.type
        rows["bias_factor"]     = lfparams.bias_factor
        rows["bias_target"]     = lfparams.bias_target
        rows["weight"]          = lfparams.weight
        rows["max_x"]           = lfparams.max_x
        return 1

tableio.I3ConverterRegistry.register(I3CorsikaWeightConverter)
