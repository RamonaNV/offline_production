from I3Tray import *
from icecube import icetray, dataclasses
from icecube.dataclasses import I3Particle
import logging

import math
from collections import namedtuple

# Copied form CORSIKA
particle_id_mass = {
       I3Particle.PPlus       : 1.00797,
       I3Particle.PMinus      : 1.00797,
       I3Particle.He3Nucleus  : 4.0026,
       I3Particle.He4Nucleus  : 4.0026,
       I3Particle.Be9Nucleus  : 9.0122,
       I3Particle.C12Nucleus  : 12.0112,
       I3Particle.N14Nucleus  : 14.0067,
       I3Particle.O16Nucleus  : 15.9994,
       I3Particle.O18Nucleus  : 15.9994,
       I3Particle.Ne20Nucleus : 20.183,
       I3Particle.Mg24Nucleus : 24.312,
       I3Particle.Al26Nucleus : 26.9815,
       I3Particle.Al27Nucleus : 26.9815,
       I3Particle.Si28Nucleus : 28.086,
       I3Particle.Si30Nucleus : 28.086,
       I3Particle.P32Nucleus  : 30.984,
       I3Particle.S35Nucleus  : 32.064,
       I3Particle.Ar39Nucleus : 39.948,
       I3Particle.K39Nucleus  : 39.102,
       I3Particle.Ca40Nucleus : 40.08,
       I3Particle.Sc44Nucleus : 44.956,
       I3Particle.Ti48Nucleus : 47.9,
       I3Particle.Ti44Nucleus : 47.9,
       I3Particle.Ti45Nucleus : 47.9,
       I3Particle.Ti46Nucleus : 47.9,
       I3Particle.Ti47Nucleus : 47.9,
       I3Particle.Ti49Nucleus : 47.9,
       I3Particle.Ti50Nucleus : 47.9,
       I3Particle.V49Nucleus  : 50.942,
       I3Particle.V48Nucleus  : 50.942,
       I3Particle.V50Nucleus  : 50.942,
       I3Particle.Mn54Nucleus : 54.938,
       I3Particle.Fe56Nucleus : 55.847,
}

particle_id_charge = {
       I3Particle.PPlus       : 1,
       I3Particle.PMinus      : 1,
       I3Particle.He3Nucleus  : 2,
       I3Particle.He4Nucleus  : 2,
       I3Particle.Be9Nucleus  : 4,
       I3Particle.C12Nucleus  : 6,
       I3Particle.N14Nucleus  : 7,
       I3Particle.O18Nucleus  : 8,
       I3Particle.O16Nucleus  : 8,
       I3Particle.Ne20Nucleus : 10 ,
       I3Particle.Mg24Nucleus : 12 ,
       I3Particle.Si28Nucleus : 14 ,
       I3Particle.Al26Nucleus : 13 ,
       I3Particle.Al27Nucleus : 13 ,
       I3Particle.Si30Nucleus : 14 ,
       I3Particle.P32Nucleus  : 15 ,
       I3Particle.S35Nucleus  : 16 ,
       I3Particle.K39Nucleus  : 19 ,
       I3Particle.Ar39Nucleus : 18 ,
       I3Particle.Ca40Nucleus : 20 ,
       I3Particle.Sc44Nucleus : 21 ,
       I3Particle.Ti48Nucleus : 22 ,
       I3Particle.Ti50Nucleus : 22 ,
       I3Particle.Ti49Nucleus : 22 ,
       I3Particle.Ti47Nucleus : 22 ,
       I3Particle.Ti45Nucleus : 22 ,
       I3Particle.Ti44Nucleus : 22 ,
       I3Particle.Ti46Nucleus : 22 ,
       I3Particle.V50Nucleus  : 23 ,
       I3Particle.V48Nucleus  : 23 ,
       I3Particle.V49Nucleus  : 23 ,
       I3Particle.Mn54Nucleus : 25 ,
       I3Particle.Fe56Nucleus : 26 ,
}



def CalcAreaSum(cylinderRadius,cylinderLength,thetaMin,thetaMax):
	 """
	  Integrated area*solid_angle for a cylindrical sampling surface
	 """
	 from icecube.dataclasses import I3Constants
	 areaSum = I3Constants.pi * I3Constants.pi * cylinderRadius * cylinderRadius * \
					 (math.cos(thetaMin)*math.cos(thetaMin) - math.cos(thetaMax)*math.cos(thetaMax)) + \
			   2. * I3Constants.pi * cylinderRadius * cylinderLength * (thetaMax - thetaMin) - \
			I3Constants.pi * cylinderRadius * cylinderLength * (math.sin(2.*thetaMax) - math.sin(2.*thetaMin))

	 areaSum *= I3Units.steradian/I3Units.meter2

	 """
	   if thetaMin = 0. and thetaMax_ = pi/2 (90 degrees) :
	   double areaSum = I3Constants::pi * I3Constants::pi * I3Units::steradian *
						(cylinderRadius_) * (cylinderRadius_ + cylinderLength_)/I3Units::meter2;
	 """
	 return areaSum


def CalcFluxSum(energyPrimaryMin,energyPrimaryMax,primarySpectralIndex,CRMass):
	 """
	  Integrated total flux
	 """
	 # Calculate FluxSum (energy-integrated primary flux) for each primary particle
	 if primarySpectralIndex == -1 :# if E^-1 use log
		   fluxSum = math.log(energyPrimaryMax/energyPrimaryMin)

	 else: # if not E^-1 integrate power law
		   fluxSum = ( 1.0/ (1. + primarySpectralIndex) ) * \
				 ( pow( energyPrimaryMax*CRMass / I3Units.TeV, 1.+primarySpectralIndex ) - \
				   pow( energyPrimaryMin*CRMass / I3Units.TeV, 1.+primarySpectralIndex) )

	 return fluxSum

def CalcHoerandelFluxSum(emin,dslope=0.0):
    """
    Calculation of FluxSum for Hoerandel spectrum
    """
    import math

    norm = [ 0.0873, 0.0571, 0.00208, 0.000474, 0.000895, 0.0106, 0.00235, 0.0157, 0.000328, 0.0046,
       0.000754, 0.00801, 0.00115, 0.00796, 0.00027, 0.00229, 0.000294, 0.000836, 0.000536, 0.00147,
       0.000304, 0.00113, 0.000631, 0.00136, 0.00135, 0.0204 ]
    nmax = len(norm)

    gamma = [ 2.71, 2.64, 2.54, 2.75, 2.95, 2.66, 2.72, 2.68, 2.69, 2.64, 2.66, 2.64, 2.66, 2.75,
       2.69, 2.55, 2.68, 2.64, 2.65, 2.7, 2.64, 2.61, 2.63, 2.67, 2.46, 2.59 ]

    crs = [ 1.00797, 4.0026, 6.939, 9.0122, 10.811, 12.0112, 14.0067, 15.9994, 18.9984, 20.183,
       22.9898, 24.312, 26.9815, 28.086, 30.984, 32.064, 35.453, 39.948, 39.102, 40.08, 44.956,
       47.9, 50.942, 51.996, 54.938, 55.847 ]

    integral = 0.0
    for i in range(nmax):
      prwght = round(crs[i])
      gamma[i] += dslope
      integral += norm[i] * math.pow( (prwght*emin/I3Units.TeV), (1-gamma[i]) ) / (gamma[i]-1)

    return integral




def FiveComponent(nshowers, emin, emax, normalization=[10., 5., 3., 2., 1.], gamma=[-2.]*5,
    LowerCutOffType="EnergyPerNucleon", UpperCutOffType="EnergyPerParticle"):
	"""
	Mimic 5-component dCORSIKA with configurations for 5 single-primary runs.
	
	:param nshowers: Total number of showers to generate
	:param emin: Minimum primary energy [GeV]
	:param emax: Maximum primary energy [GeV]
	:param normalization: Ratios of the total number of generated showers in each component
	:param gamma: Spectral index of each component
	:param LowerCutOffType: If EnergyPerNucleon, make lower bound of energy range proportional to A
	:param UpperCutOffType: If EnergyPerNucleon, make upper bound of energy range proportional to A
	"""
	Component = namedtuple('PowerLawComponent', ['NEvents', 'PrimaryType', 'EnergyPrimaryMin', 'EnergyPrimaryMax', 'PrimarySpectralIndex'])
	
	def integrate(emin, emax, index):
		if index == -1:
			return math.log(emax/emin)
		elif index < -1:
			g = index + 1
			return (emax**g - emin**g)/g
	
	ptype = [getattr(I3Particle, p) for p in ('PPlus', 'He4Nucleus', 'N14Nucleus', 'Al27Nucleus', 'Fe56Nucleus')]
	corsika_ptype = [14, 402, 1407, 2713, 5626]
	if LowerCutOffType == "EnergyPerNucleon":
		emin = [emin*round(particle_id_mass[pt]) for pt in ptype]
	else:
		emin = [emin]*len(ptype)
	if UpperCutOffType == "EnergyPerNucleon":
		emax = [emax*round(particle_id_mass[pt]) for pt in ptype]
	else:
		emax = [emax]*len(ptype)
	integral = [n*integrate(elo/I3Units.TeV, ehi/I3Units.TeV, index) for n, index, elo, ehi in zip(normalization, gamma, emin, emax)]
	total = sum(integral)
	normalization = [part/total for n, part in zip(normalization, integral)]
	
	
	return [Component(max((int(round(n*nshowers)), 1)), pt, elo, ehi, index) for n, pt, index, elo, ehi in zip(normalization, corsika_ptype, gamma, emin, emax)]

class Corsika5CompWeightModule(icetray.I3Module):

     def __init__(self,context):
         icetray.I3Module.__init__(self, context)
         self.logger = logging.getLogger(self.__class__.__name__)
         self.AddOutBox("OutBox");

         self.AddParameter("NEvents", "Number of generated showers in CORSIKA", -1) # no sane default
         self.AddParameter("Name", "Name to give the frame", "CorsikaWeightMap")
         self.AddParameter("ThetaMin", "Primary cosmic ray minimum zenith angle", 0.0)
         self.AddParameter("ThetaMax", "Primary cosmic ray maximum zenith angle", 89.9*I3Units.deg)
         self.AddParameter("MCTreeName", "Name of MCTree List to get from frame", "I3MCTree")
         self.AddParameter("CylinderRadius", "Radius of Generation Cylinder", -1)
         self.AddParameter("CylinderLength", "Height of Generation Cylinder", -1)
         self.AddParameter("EnergyPrimaryMin", "Primary minimum energy (from CORSIKA steering file)", 0.0)
         self.AddParameter("EnergyPrimaryMax", "Primary maximum energy (from CORSIKA steering file)", 0.0)

         self.AddParameter("PrimaryNormalizationH", 
                "Primary normalization of generated nuclei (i.e. steering parameter PNORM)", 1.0)
         self.AddParameter("PrimaryNormalizationHe", 
                "Primary normalization of generated nuclei (i.e. steering parameter PNORM)", 1.0)
         self.AddParameter("PrimaryNormalizationCNO", 
                "Primary normalization of generated nuclei (i.e. steering parameter PNORM)", 1.0)
         self.AddParameter("PrimaryNormalizationMgAlSi", 
                "Primary normalization of generated nuclei (i.e. steering parameter PNORM)", 1.0)
         self.AddParameter("PrimaryNormalizationFe", 
                "Primary normalization of generated nuclei (i.e. steering parameter PNORM)", 1.0)

         self.AddParameter("PrimarySpectralIndexH", 
                "Primary spectral index of generated nuclei (i.e. steering parameter PGAM)", -2.0)
         self.AddParameter("PrimarySpectralIndexHe", 
                "Primary spectral index of generated nuclei (i.e. steering parameter PGAM)", -2.0)
         self.AddParameter("PrimarySpectralIndexCNO", 
                "Primary spectral index of generated nuclei (i.e. steering parameter PGAM)", -2.0)
         self.AddParameter("PrimarySpectralIndexMgAlSi", 
                "Primary spectral index of generated nuclei (i.e. steering parameter PGAM)", -2.0)
         self.AddParameter("PrimarySpectralIndexFe", 
                "Primary spectral index of generated nuclei (i.e. steering parameter PGAM)", -2.0)

         self.AddParameter("Spric", "Separate primary energy cutoff (i.e. steering parameter SPRIC)", True)
         self.AddParameter("OverSamplingFactor", "Particle oversampling factor", 1.0)

     def Configure(self):
         self.corsikaShowers   =  self.GetParameter("NEvents")
         assert self.corsikaShowers >= 0, "You must specify NEvents"

         self.weight_name      =  self.GetParameter("Name")
         self.thetaMin         =  self.GetParameter("ThetaMin") / I3Units.rad
         self.thetaMax         =  self.GetParameter("ThetaMax") / I3Units.rad
         self.mcTreeName       =  self.GetParameter("MCTreeName")
         self.cylinderRadius   =  self.GetParameter("CylinderRadius") * I3Units.meter
         self.cylinderLength   =  self.GetParameter("CylinderLength") * I3Units.meter
         self.energyPrimaryMin =  self.GetParameter("EnergyPrimaryMin") * I3Units.GeV
         self.energyPrimaryMax =  self.GetParameter("EnergyPrimaryMax") * I3Units.GeV

         self.primaryNormalization = []
         self.primaryNormalization.append( self.GetParameter("PrimaryNormalizationH") )
         self.primaryNormalization.append( self.GetParameter("PrimaryNormalizationHe") )
         self.primaryNormalization.append( self.GetParameter("PrimaryNormalizationCNO") )
         self.primaryNormalization.append( self.GetParameter("PrimaryNormalizationMgAlSi") )
         self.primaryNormalization.append( self.GetParameter("PrimaryNormalizationFe") )

         self.primarySpectralIndex = []
         self.primarySpectralIndex.append( self.GetParameter("PrimarySpectralIndexH") )
         self.primarySpectralIndex.append( self.GetParameter("PrimarySpectralIndexHe") )
         self.primarySpectralIndex.append( self.GetParameter("PrimarySpectralIndexCNO") )
         self.primarySpectralIndex.append( self.GetParameter("PrimarySpectralIndexMgAlSi") )
         self.primarySpectralIndex.append( self.GetParameter("PrimarySpectralIndexFe") )

         self.spric          =  self.GetParameter("Spric")
         self.oversampling   =  self.GetParameter("OverSamplingFactor") 

         self.areaSum = CalcAreaSum(self.cylinderRadius,self.cylinderLength,self.thetaMin,self.thetaMax)
         self.fluxSum = [0.]*5;

         self.particle_ids = [
                I3Particle.PPlus,       # H
                I3Particle.He4Nucleus,   # He
                I3Particle.N14Nucleus,  # N
                I3Particle.Al27Nucleus,  # Al
                I3Particle.Fe56Nucleus   # Fe 
                ]

         # normalize primary normalization to 1
         normsum = sum(self.primaryNormalization)
         self.primaryNormalization = map(lambda x:x/normsum, self.primaryNormalization)

         # loop through the 5 mass group
         normsum = 0.; 

         # Calculate FluxSum (energy-integrated primary flux) for each primary particle
         for i in range(5):
            # cosmic ray mass number
            # if SPRIC == 1 then energy/nucleon is assumed
            # if SPRIC != 1 the energy/particle is assumed
            if self.spric: 
               cr_mass = round(particle_id_mass[self.particle_ids[i]])
            else:
               cr_mass = 1.0

            self.fluxSum[i] = CalcFluxSum(
                                      self.energyPrimaryMin,self.energyPrimaryMax,
                                      self.primarySpectralIndex[i],cr_mass)
            normsum +=  self.fluxSum[i]*self.primaryNormalization[i] 

         assert (normsum > 0.), "normalization factor not properly initialized"
         self.events = [0.]*5;
         for i in range(5):
            # normalize primary normalization factor
            self.primaryNormalization[i] *= self.fluxSum[i]/normsum;

            # calculate scaled numbers of showers from a given primary particle
            self.events[i] = round(self.corsikaShowers * self.primaryNormalization[i]);


     def DAQ(self,frame):
         """
         Process the frames
         """
         # Get the I3MCTree from the frame
         if not self.mcTreeName in frame:
            self.logger.error( "Frame does not contain MCTree: %s",self.mcTreeName)
         mctree = frame[self.mcTreeName]

         # Get list of primaries
         primaries = mctree.primaries
         if len(primaries) > 1:
            self.logger.warn( "%s contains more than one primary. Will attempt to determine from weights",self.mcTreeName)
            if self.weight_name in frame: # Try to determine primary from weights
               weightdict = frame[self.weight_name]
               primary = findPrimary(primaries,weightdict)

         elif len(primaries) == 0:
             self.logger.error( "%s contains no primary",self.mcTreeName)
         else:
            primary = primaries[0]
         # Get the primary particle energy: i.e. the first of the list
         cosmicParticleEnergy = primary.energy/I3Units.GeV

         # Get the primary particle type: i.e. the first of the list
         particle_type         =  primary.type
         particle_index        =  self.particle_ids.index(particle_type)
         primarySpectralIndex  =  self.primarySpectralIndex[particle_index]
         if self.spric: 
               cr_mass = particle_id_mass[particle_type]
         else:
               cr_mass = 1.0

         # Instantiate a Map where to store weighting info
         if self.weight_name in frame: # If weightmap is there, set FluxSum and TimeScale
            weightdict           = frame[self.weight_name]
            nevents              = weightdict["NEvents"]
            if "AreaSum" not in weightdict:
            	weightdict["AreaSum"] = CalcAreaSum(self.cylinderRadius,self.cylinderLength,self.thetaMin,self.thetaMax)
            areaSum              = weightdict["AreaSum"]  
            energyPrimaryMin     = weightdict["EnergyPrimaryMin"]     
            energyPrimaryMax     = weightdict["EnergyPrimaryMax"]     
            oldPrimarySpectralIndex = weightdict["PrimarySpectralIndex"] 

            primarySpectralIndex               = self.primarySpectralIndex[particle_index]
            weightdict["PrimarySpectralIndex"] = primarySpectralIndex

            # This CRMass is explicitly set to 1.0 since the energy range already takes this into account
            fluxSum              = CalcFluxSum(energyPrimaryMin,energyPrimaryMax,
                                               primarySpectralIndex,1.0)

            weightdict["FluxSum"]   = fluxSum
            weightdict["TimeScale"] = nevents / (areaSum * fluxSum)
            weightdict["OverSampling"] = self.oversampling 

         else:
            weightdict                         = dataclasses.I3MapStringDouble()
            particle_index                     = self.particle_ids.index(particle_type)
            primarySpectralIndex               = self.primarySpectralIndex[particle_index]

            # Start filling the Map
            weightdict["PrimaryType"]          = primary.type
            weightdict["PrimarySpectralIndex"] = primarySpectralIndex
            
            weightdict["CylinderRadius"]       = self.cylinderRadius
            weightdict["CylinderLength"]       = self.cylinderLength
            weightdict["AreaSum"]              = self.areaSum
            weightdict["EnergyPrimaryMin"]     = self.energyPrimaryMin
            weightdict["EnergyPrimaryMax"]     = self.energyPrimaryMax

            weightdict["NEvents"]              = self.events[particle_index]
            fluxSum                            = CalcFluxSum(self.energyPrimaryMin,self.energyPrimaryMax,
                                                               primarySpectralIndex,cr_mass)
            weightdict["FluxSum"]              = fluxSum
            weightdict["TimeScale"]            = self.events[particle_index]  / (self.areaSum * fluxSum)
            weightdict["OverSampling"]         = self.oversampling 



         # Calculate weight to flat spectrum for each primary particle
         energyFactor = pow( cosmicParticleEnergy / I3Units.TeV , primarySpectralIndex )
         oldEnergyFactor = pow( cosmicParticleEnergy / I3Units.GeV , primarySpectralIndex )

         # CorsikaWeight has units of (GeV m^2 sec sr)
         CorsikaWeight           = 1.0 / energyFactor
         weightdict["Weight"]    = CorsikaWeight

         oldCorsikaWeight        = 1.0 / oldEnergyFactor
         weightdict["OldWeight"] = oldCorsikaWeight

         weightdict["PrimaryEnergy"] = cosmicParticleEnergy/I3Units.GeV

         if "PrimaryType" not in weightdict:
            weightdict["PrimaryType"]  = primary.type

         if self.weight_name in frame: # Remove old weightmap
            del frame[self.weight_name]
         frame[self.weight_name] = weightdict
         if not frame.Has("PrimaryParticle"):
            frame["PrimaryParticle"] = primary

         self.PushFrame(frame,"OutBox")




class CorsikaWeightModule(icetray.I3Module):
     """
     Sets TimeScale and FluxSum for unweighted CORSIKA
     """
     def __init__(self,context):
         icetray.I3Module.__init__(self, context)
         self.logger = logging.getLogger(self.__class__.__name__)
         self.AddOutBox("OutBox");

         self.AddParameter("NEvents", "Number of generated showers in CORSIKA", -1) # no sane default
         self.AddParameter("Name", "Name to give the frame", "CorsikaWeightMap")
         self.AddParameter("ThetaMin", "Primary cosmic ray minimum zenith angle", 0.0)
         self.AddParameter("ThetaMax", "Primary cosmic ray maximum zenith angle", 89.9*I3Units.deg)
         self.AddParameter("MCTreeName", "Name of MCTree List to get from frame", "I3MCTree")
         self.AddParameter("Spric", "Separate primary energy cutoff (i.e. steering parameter SPRIC)", True)
         self.AddParameter("CylinderRadius", "Radius of Generation Cylinder", -1)
         self.AddParameter("CylinderLength", "Height of Generation Cylinder", -1)
         self.AddParameter("OverSamplingFactor", "Oversampling factor", 1)
 

     def Configure(self):
         self.mcTreeName       =  self.GetParameter("MCTreeName")
         self.coriskaShowers   =  self.GetParameter("NEvents")
         self.weight_name      =  self.GetParameter("Name")
         self.spric            = self.GetParameter("Spric")
         self.thetaMin         = self.GetParameter("ThetaMin")
         self.thetaMax         = self.GetParameter("ThetaMax")
         self.cylinderRadius   = self.GetParameter("CylinderRadius")
         self.cylinderLength   = self.GetParameter("CylinderLength")
         self.oversampling     = self.GetParameter("OverSamplingFactor")



     def DAQ(self,frame):
         """
         Process the frames
         """
         # Get the I3MCTree from the frame
         if self.mcTreeName not in frame:
             self.logger.error( "Frame does not contain %s",self.mcTreeName)
             self.PushFrame(frame,"OutBox")
             return
         mctree = frame[self.mcTreeName]

         # Get list of primaries
         primaries = mctree.primaries
         if len(primaries) > 1:
             self.logger.warning( "!!!!%s contains more than one primary",self.mcTreeName)
             self.PushFrame(frame,"OutBox")
             return
         elif len(primaries) == 0:
             self.logger.error( "!!!!%s contains no primary",self.mcTreeName)
             self.PushFrame(frame,"OutBox")
             return

         # Get the primary particle energy: i.e. the first of the list
         cosmicParticleEnergy = primaries[0].energy

         # Get the primary particle type: i.e. the first of the list
         particle_type        =  primaries[0].type
         cr_mass              =  1.0
         
         # Instantiate a Map where to store weighting info
         if self.weight_name in frame: # If weightmap is ther set FlusSum and TimeScale
            weightdict           = frame[self.weight_name]
            nevents              = weightdict["NEvents"]
            if "AreaSum" not in weightdict:
            	weightdict["AreaSum"] = CalcAreaSum(self.cylinderRadius,self.cylinderLength,self.thetaMin,self.thetaMax)
            areaSum              = weightdict["AreaSum"]  
            energyPrimaryMin     = weightdict["EnergyPrimaryMin"]     
            energyPrimaryMax     = weightdict["EnergyPrimaryMax"]     
            primarySpectralIndex = weightdict["PrimarySpectralIndex"] 
            fluxSum              = CalcHoerandelFluxSum(energyPrimaryMin)

            weightdict["FluxSum"]       = fluxSum
            weightdict["TimeScale"]     = nevents / (areaSum * fluxSum)
            weightdict["OverSampling"]  = self.oversampling 
            weightdict["Weight"]        = 1.0
            weightdict["PrimaryEnergy"] = cosmicParticleEnergy/I3Units.GeV
            weightdict["ParticleType"]  = particle_type
            weightdict["Polygonato"]    = 1.0


            if self.weight_name in frame: # Remove o<D->>ld weightmap
                del frame[self.weight_name]
            frame[self.weight_name] = weightdict

            self.PushFrame(frame,"OutBox")


class PolygonatoWeightModule(icetray.I3Module):
  
  def __init__(self,ctx):
      icetray.I3Module.__init__(self,ctx)
      self.logger = logging.getLogger(self.__class__.__name__)
      self.AddParameter("WeightMapName","Name of object where weights are stored","CorsikaWeightMap")
      self.AddParameter("MCTreeName","Name of I3MCTree object in frame","I3MCTree")
      self.AddOutBox("OutBox")

  def Configure(self):
      self.weightmapname = self.GetParameter("WeightMapName")
      self.mctreename    = self.GetParameter("MCTreeName")
      
  def DAQ(self,frame):
      mctree    = frame[self.mctreename]
      weightmap = frame[self.weightmapname]
      primaries = mctree.primaries
      primary   = primaries[0]
      if len(mctree.primaries) > 1:
          self.logger.warn( "%s contains more than one primary. Will attempt to determine from weights",self.mctreename)
          primary = findPrimary(mctree.primaries,weightmap)

      tFlux    = 0.0
      pslope   = 0.0
      nevents  = 0.0
      zPrimary = particle_id_charge[primary.type]
      ePrimary = primary.energy/I3Units.TeV

      if primary.type == I3Particle.PPlus:
          tFlux = 8.73e-2
          pslope = pow(ePrimary,-2.71)

      elif primary.type == I3Particle.He4Nucleus:
          tFlux = 5.71e-2
          pslope = pow(ePrimary,-2.64)

      elif primary.type == I3Particle.N14Nucleus:
          tFlux = 3.24e-2
          pslope = pow(ePrimary,-2.67)

      elif primary.type == I3Particle.Al27Nucleus:
          tFlux = 3.16e-2
          pslope = pow(ePrimary,-2.65)

      elif primary.type == I3Particle.Fe56Nucleus:
          tFlux = 2.18e-2
          pslope = pow(ePrimary,-2.575)
      else:
          tFlux = 0.0;
          self.logger.error( "Error: invalid primary type %s", primary.type)

      rigidity   = pow(1.+ pow(ePrimary/(4.49e3*zPrimary),1.9),-2.1/1.9)
      polygonato = tFlux*pslope*rigidity

      weightmap["Polygonato"] = polygonato
      del frame[self.weightmapname]
      frame.Put(self.weightmapname,weightmap)

      self.PushFrame(frame)


def findPrimary(primaries,weightdict):
    """
    Sometimes you have to reverse engineer the weights to extract the primary from the frame
    """
    oldPrimarySpectralIndex = weightdict["PrimarySpectralIndex"] 
    eFactor = lambda x: 1.0/pow( x / I3Units.TeV , oldPrimarySpectralIndex)
    for p in primaries:
        if abs(eFactor(p.energy) - weightdict["Weight"]) < 1.0e-6: 
           if "PrimaryType" in weightdict.keys():
              if abs(float(p.type) - weightdict["PrimaryType"]) < 1.0e-6:
                 return p
           else: return p
    return None

