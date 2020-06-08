#!/usr/bin/env python
import os
import subprocess
import threading
import time
import logging

from os.path import expandvars
from icecube import icetray, dataclasses

import math

def taskset(pid,tt=None):
    # get/set the taskset affinity for pid
    # uses a binary number string for the core affinity
    l = ['/bin/taskset','-p']
    if tt:
        l.append(hex(int(tt,2))[2:])
    l.append(str(pid))
    p = subprocess.Popen(l,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output = p.communicate()[0].split(':')[-1].strip()
    if not tt:
        return bin(int(output,16))[2:]

def tasksetInUse():
    # check for cpu affinity (taskset)
    try:
        num_cpus = reduce(lambda b,a: b+int('processor' in a),open('/proc/cpuinfo').readlines(),0)
        affinity = taskset(os.getpid())
        if len(affinity) < num_cpus:
            return True
        for x in affinity[:num_cpus]:
            if x != '1':
                return True
        return False
    except Exception:
        return False

def resetTasksetThreads(main_pid):
    # reset thread taskset affinity
    time.sleep(60)
    try:
        num_cpus = reduce(lambda b,a: b+int('processor' in a),open('/proc/cpuinfo').readlines(),0)
        tt = '1'*num_cpus
        p = subprocess.Popen(['/bin/ps','-Lo','tid','--no-headers','%d'%main_pid],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        for tid in p.communicate()[0].split():
            tid = tid.strip()
            if tid:
                taskset(tid,tt)
    except Exception:
        pass

def LoadCascadeTables(IceModel = "SpiceMie", TablePath = "/data/sim/sim-new/spline-tables"):
    from icecube import photonics_service
    if IceModel == "Spice1":
        amplitudetable = TablePath+'/ems_spice1_z20_a10.abs.fits'
        timingtable = TablePath+'/ems_spice1_z20_a10.prob.fits'
    elif IceModel == "SpiceMie":
        amplitudetable = TablePath+'/ems_mie_z20_a10.abs.fits'
        timingtable = TablePath+'/ems_mie_z20_a10.prob.fits'
    elif IceModel == "SpiceMieNoHoleIce":
        amplitudetable = TablePath+'/NoHoleIceCascades_250_z20_a10.abs.fits'
        timingtable = TablePath+'/NoHoleIceCascades_250_z20_a10.prob.fits'
    else:
        raise RuntimeError("Unknown ice model: %s", IceModel)
    
    logging.debug("Loading cascade tables : ")
    logging.debug("  amp  = ", amplitudetable)
    logging.debug("  time = ", timingtable)

    cascade_service = photonics_service.I3PhotoSplineService(
        amplitudetable=amplitudetable,
        timingtable=timingtable,
        timingSigma=0.,
        maxRadius=800.*icetray.I3Units.meter)

    return cascade_service

@icetray.traysegment
def PropagatePhotons(tray, name,
    GCDFile,
    If=lambda f:True ,
    RandomService = None,
    KeepIndividualMaps = False,
    HybridMode = False,
    IgnoreMuons = False,
    IgnoreCascades = False,
    UseGPUs = False,
    UseAllCPUCores = False,
    KeepSlicedMCTree = False,
    IceModel = "spice_3.2",
    CascadeService = None,
    IceModelLocation = None,
    UseCascadeExtension = True,
    UseGeant4=False, 
    CrossoverEnergyEM=None, 
    CrossoverEnergyHadron=None, 
    UnshadowedFraction = 1.0, #changed 2014-10-16 to IC86 nominal preset, IC79 used 0.9
    DOMOversizeFactor=5.0,
    HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"),
    InputMCTree="I3MCTree",
    UseI3PropagatorService = False,
    OutputPESeriesMapName="I3MCPESeriesMap",
    OutputPhotonSeriesName=None,
):
    """ This traysegment offers multiple tweaks for adapted processing in different energy ranges,
    for example GEANT4 in conjunction with Parametrizations for the treatment for lowest energies 
    and a HybridMode with the use of tables for the treatment of high energies.
    In any case, please refer to the documentation of clsim to find suitable settings for your 
    simulation needs """

    from I3Tray import I3Units
    from icecube import icetray, dataclasses, dataio
    from icecube import phys_services, sim_services

    from icecube import clsim, photonics_service
    from os import listdir
    from os.path import isdir

    
    if IgnoreMuons and not HybridMode:
        raise RuntimeError("Can currently only ignore muons in hybrid mode")

    clsimIceModel = None
    if IceModelLocation is None:
       IceModelLocation = expandvars("$I3_BUILD/ice-models/resources/models")
    if isinstance(IceModel, clsim.I3CLSimMediumProperties):
        if HybridMode:
            raise RuntimeError("Cannot use custom ice models in hybrid mode")
        clsimIceModel = IceModel
    elif IceModel == "Spice1":
        clsimIceModel = expandvars(IceModelLocation+"/spice_1")
    elif IceModel == "SpiceMie":
        clsimIceModel = expandvars(IceModelLocation+"/spice_mie")
    elif IceModel == "SpiceLea":
        clsimIceModel = expandvars(IceModelLocation+"/spice_lea")
    else:
        for d in listdir(IceModelLocation):
           if isdir(expandvars(IceModelLocation+"/"+d)) and IceModel.lower() == d.lower():
              clsimIceModel = expandvars(IceModelLocation+"/"+d)
              break
        if not clsimIceModel:
           raise RuntimeError("Unknown ice model: %s", IceModel)

    if HybridMode and IceModel not in ("Spice1","SpiceMie"):
            raise RuntimeError("Can only use Spice1 and SpiceMie in hybrid mode. photon tables do not support ice anisotropy at this time.")

    if (not IgnoreCascades) and HybridMode:
        if CascadeService is None:
            logging.warning("*** no cascades tables provided. Loading tables for", IceModel)
            
            # If we can see CVMFS, we'll get the splines from there.
            #  Note : when available, switch icecube.wisc.edu for icecube.opensciencegrid.org
            UseSplinesFromCVMFS = os.path.isdir("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines")

            if(UseSplinesFromCVMFS):
                TablePath="/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines"
            else:
                TablePath="/data/sim/sim-new/spline-tables"                
            
            logging.info("Using splines from CVMFS: ", UseSplinesFromCVMFS)

            # Work out which splines to use based on ice model preferences
            if(HoleIceParameterization==expandvars('$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm')):
                CascadeModel=IceModel
            elif(HoleIceParameterization==expandvars('$I3_SRC/ice-models/resources/models/angsens/as.nominal')):
                if IceModel=="SpiceMie":
                    CascadeModel="SpiceMieNoHoleIce"
                else:
                    raise RuntimeError("No no-hole-ice spline for %s", IceModel)
            else: 
                raise RuntimeError("No spline for %s with hole ice param %s", IceModel, HoleIceParameterization)

            cascade_service = LoadCascadeTables(IceModel=CascadeModel, TablePath=TablePath)
        else:
            cascade_service = CascadeService
    else:
        cascade_service = None

    if HybridMode:
        if OutputPhotonSeriesName is not None:
            raise RuntimeError("saving photons is not supported in hybrid mode")
        if UseGeant4:
            raise RuntimeError("Geant4 not supported in hybrid mode")
        if ((CrossoverEnergyEM is not None) or (CrossoverEnergyHadron is not None)):
            raise RuntimeError("CrossoverEnergyEM or CrossoverEnergyHadron not supported in hybrid mode")

        # split the MCTree into a cascade-only and a track-only version
        tray.AddModule("I3MCTreeHybridSimulationSplitter", name+"_splitMCTree",
            InputMCTreeName=InputMCTree,
            OutputMCTreeNameTracks=InputMCTree+"Tracks",
            OutputMCTreeNameCascades=InputMCTree+"Cascades")
        
        tray.AddModule("I3TauSanitizer", name+"_sanitize_taus",
                       InputMCTreeName = InputMCTree+"Tracks",
                       OutputMCTreeName = InputMCTree+"Tracks") # overwrite the input

        if not IgnoreMuons:
            if UseGPUs:
                DoNotParallelize=False
            else:
                DoNotParallelize=not UseAllCPUCores
                threading.Thread(target=resetTasksetThreads,args=(os.getpid(),)).start()
            logging.debug('tasksetInUse = ',tasksetInUse())
            logging.debug('DoNotParallelize = ',DoNotParallelize)

            # simulate tracks (with clsim)
            tray.AddSegment(clsim.I3CLSimMakeHits, name+"_makeCLSimHits",
                PhotonSeriesName = None,
                MCTreeName = InputMCTree+"Tracks",
                OutputMCTreeName = InputMCTree+"Tracks_sliced",
                MCPESeriesName = OutputPESeriesMapName + "Tracks",
                UseI3PropagatorService = UseI3PropagatorService,
                RandomService = RandomService,
                UnshadowedFraction=UnshadowedFraction,
                DoNotParallelize=DoNotParallelize,
                UseGeant4=False, # never use this with Geant4!
                UseGPUs=UseGPUs,
                UseCPUs=not UseGPUs,
                IceModelLocation=clsimIceModel,
                DOMOversizeFactor=DOMOversizeFactor,
                UseCascadeExtension=UseCascadeExtension,
                GCDFile=GCDFile,
                DisableTilt=True)

            tray.AddModule("Delete", name+"_cleanup_clsim_sliced_MCTree",
                Keys = [InputMCTree+"Tracks_sliced"])

        if not IgnoreCascades:
            tray.AddModule("I3PhotonicsHitMaker", name+"_hitsFromTheTable",
                CascadeService = cascade_service,
                TrackService = None, # tracks are handled by clsim
                UnshadowedFraction = UnshadowedFraction,
                Input = InputMCTree+"Cascades",
                Output = OutputPESeriesMapName + "Cascades",
                RandomService = RandomService
                )

        MCPEsToCombine = []
        if not IgnoreMuons:
            MCPEsToCombine.append(OutputPESeriesMapName + "Tracks")
        if not IgnoreCascades:
            MCPEsToCombine.append(OutputPESeriesMapName + "Cascades")

        # combine the resulting I3MCHitSeriesMaps
        tray.AddModule("I3CombineMCPE", name+"_combine_pes",
            InputResponses = MCPEsToCombine,
            OutputResponse = OutputPESeriesMapName)

        if not KeepIndividualMaps:
            # delete the original maps and the split I3MCTrees
            tray.AddModule("Delete", name+"_cleanup_peseriesmaps",
                Keys = MCPEsToCombine)

            tray.AddModule("Delete", name+"_cleanup_MCTree",
                Keys=[InputMCTree+"Tracks", InputMCTree+"Cascades"])

    else:
        # non-hybrid clsim-only simulation
        # If we're using Geant4, we do NOT want the taus to be dark.
        if not UseGeant4:
            tray.AddModule("I3TauSanitizer", name+"_sanitize_taus",
                InputMCTreeName = InputMCTree,
                OutputMCTreeName = InputMCTree) # overwrite the input

        if UseGPUs:
            DoNotParallelize=False
        else:
            DoNotParallelize=not UseAllCPUCores
            threading.Thread(target=resetTasksetThreads,args=(os.getpid(),)).start()
        logging.debug('tasksetInUse = %s' % tasksetInUse())
        logging.debug('DoNotParallelize = %s' % DoNotParallelize)

        # simulate tracks (with clsim)
        tray.AddSegment(clsim.I3CLSimMakeHits, name+"_makeCLSimHits",
            PhotonSeriesName = OutputPhotonSeriesName,
            MCTreeName = InputMCTree,
            MCPESeriesName = OutputPESeriesMapName,
            UseI3PropagatorService = UseI3PropagatorService,
            RandomService = RandomService,
            UnshadowedFraction = UnshadowedFraction,
            DoNotParallelize = DoNotParallelize,
            UseGeant4=UseGeant4,
            CrossoverEnergyEM=CrossoverEnergyEM,
            CrossoverEnergyHadron=CrossoverEnergyHadron,
            UseGPUs=UseGPUs,
            UseCPUs=not UseGPUs,
            DOMOversizeFactor=DOMOversizeFactor,
            HoleIceParameterization=HoleIceParameterization,
            IceModelLocation=clsimIceModel,
            GCDFile=GCDFile,
            UseCascadeExtension=UseCascadeExtension)



