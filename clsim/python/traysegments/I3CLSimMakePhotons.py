#
# Copyright (c) 2011, 2012
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
# 
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 
# 
# $Id: I3CLSimMakePhotons.py 140605 2015-12-18 20:19:05Z claudio.kopper $
# 
# @file I3CLSimMakePhotons.py
# @version $Revision: 140605 $
# @date $Date: 2015-12-18 15:19:05 -0500 (Fri, 18 Dec 2015) $
# @author Claudio Kopper
#

from __future__ import print_function

import string
import numpy

from os.path import expandvars, exists, isdir, isfile
import tempfile

from icecube import icetray, dataclasses, dataio
from icecube.icetray import logging

from icecube.clsim import GetDefaultParameterizationList
from icecube.clsim import GetFlasherParameterizationList
from icecube.clsim import GetHybridParameterizationList

from .common import setupPropagators, setupDetector, configureOpenCLDevices

# use this instead of a simple "@icetray.traysegment" to support
# ancient versions of IceTray that do not have tray segments.
def unchanged(func): return func
my_traysegment = icetray.traysegment if hasattr(icetray, "traysegment") else unchanged
@my_traysegment
def I3CLSimMakePhotons(tray, name,
                       GCDFile,
                       UseCPUs=False,
                       UseGPUs=True,
                       UseOnlyDeviceNumber=None,
                       UseCUDA=False,
                       MCTreeName="I3MCTree",
                       OutputMCTreeName=None,
                       FlasherInfoVectName=None,
                       FlasherPulseSeriesName=None,
                       PhotonSeriesName="PhotonSeriesMap",
                       MCPESeriesName="MCPESeriesMap",
                       RandomService=None,
                       IceModelLocation=expandvars("$I3_BUILD/ice-models/resources/models/spice_mie"),
                       DisableTilt=False,
                       UnWeightedPhotons=False,
                       UnWeightedPhotonsScalingFactor=None,
                       UseI3PropagatorService=True,
                       UseGeant4=False,
                       ParticleHistory=False,
                       ParticleHistoryGranularity=20*icetray.I3Units.m,
                       CrossoverEnergyEM=None,
                       CrossoverEnergyHadron=None,
                       UseCascadeExtension=True,
                       StopDetectedPhotons=True,
                       PhotonHistoryEntries=0,
                       DoNotParallelize=False,
                       EnableDoubleBuffering=False,
                       DoublePrecision=False,
                       DOMOversizeFactor=5.,
                       UnshadowedFraction=0.9,
                       HoleIceParameterization=expandvars("$I3_BUILD/ice-models/resources/models/angsens/as.h2-50cm"),
                       WavelengthAcceptance=None,
                       DOMRadius=0.16510*icetray.I3Units.m, # 13" diameter
                       CableOrientation=None,
                       OverrideApproximateNumberOfWorkItems=None,
                       IgnoreSubdetectors=['IceTop'],
                       ExtraArgumentsToI3CLSimClientModule=dict(),
                       If=lambda f: True
                       ):
    """Do standard clsim processing up to the I3Photon level.
    These photons still need to be converted to I3MCPEs to be usable
    for further steps in the standard IceCube MC processing chain.
    Reads its particles from an I3MCTree and writes an I3PhotonSeriesMap.

    All available OpenCL GPUs (and optionally CPUs) will
    be used. This will take over your entire machine,
    so make sure to configure your batch jobs correctly
    when using this on a cluster.
    When using nVidia cards, you can set the
    CUDA_VISIBLE_DEVICES environment variable to
    limit GPU visibility. A setting of
    CUDA_VISIBLE_DEVICES="0,3" would only use cards
    #0 and #3 and ignore cards #1 and #2. In case you are
    using a batch system, chances are this variable is already
    set. Unfortunately, there is no corresponding setting
    for the AMD driver.

    This segment assumes that MMC has been applied to the
    I3MCTree and that MMC was *NOT* run using the "-recc" option.

    :param UseCPUs:
        Turn this on to also use CPU-based devices.
        (This may potentially slow down photon generation, which
        is also done on the CPU in parallel.)
    :param UseGPUs:
        Turn this off to not use GPU-based devices.
        This may be useful if your GPU is used for display
        purposes and you don't want it to slow down.
    :param UseOnlyDeviceNumber:
        Use only a single device number, even if there is more than
        one device found matching the required description. The numbering
        starts at 0.
    :param UseCUDA:
        Use CUDA implementation of photon propagator instead of OpenCL.
    :param MCTreeName:
        The name of the I3MCTree containing the particles to propagate.
    :param OutputMCTreeName:
        A copy of the (possibly sliced) MCTree will be stored as this name.
    :param FlasherInfoVectName:
        Set this to the name of I3FlasherInfoVect objects in the frame to
        enable flasher simulation. The module will read I3FlasherInfoVect objects
        and generate photons according to assumed parameterizations.
    :param FlasherPulseSeriesName:
        Set this to the name of an I3CLSimFlasherPulseSeries object in the frame to
        enable flasher/Standard Candle simulation.
        This cannot be used at the same time as FlasherInfoVectName.
        (I3CLSimFlasherPulseSeries objects are clsim's internal flasher
        representation, if "FlasherInfoVectName" is used, the I3FlasherInfoVect
        objects are converted to I3CLSimFlasherPulseSeries objects.)
    :param PhotonSeriesName:
        Configure this to enable writing an I3PhotonSeriesMap containing
        all photons that reached the DOM surface.
    :param RandomService:
        Set this to an instance of a I3RandomService. Alternatively,
        you can specify the name of a configured I3RandomServiceFactory
        added to I3Tray using tray.AddService(). If you don't configure
        this, the default I3RandomServiceFactory will be used.
    :param IceModelLocation:
        Set this either to a directory containing a PPC-compatible
        ice description (icemodel.dat, icemodel.par and cfg.txt) or
        to a photonics ice table file. PPC-compatible ice files should
        generally lead to faster execution times on GPUs since it involves
        less interpolation between table entries (the PPC ice-specification
        is parametric w.r.t. wavelength, whereas the photonics specification
        is not).
    :param DisableTilt:
        Do not simulate ice tilt, even if the ice model directory
        provides tilt information. (Photonics-based models will never
        have tilt.)
    :param UnWeightedPhotons:
        Enabling this setting will disable all optimizations. These
        are currently a DOM oversize factor of 5 (with the appropriate
        timing correction) and a biased initial photon spectrum that
        includes the DOM spectral acceptance. Enabling this setting
        essentially means that all photons that would be generated
        in the real detector *will* actually be generated. This will siginificatly
        slow down the simulation, but the optional ``PhotonSeries``
        will contain an unweighted sample of photons that arrive
        at your DOMs. This can be useful for DOM acceptance studies.
    :param UnWeightedPhotonsScalingFactor:
        If UnWeightedPhotons is turned on, this can be used to scale
        down the overall number of photons generated. This should normally
        not be touched (it can be used when generating photon paths
        for animation). Valid range is a float >0. and <=1.
    :param StopDetectedPhotons:
        Configures behaviour for photons that hit a DOM. If this is true (the default)
        photons will be stopped once they hit a DOM. If this is false, they continue to
        propagate.
    :param PhotonHistoryEntries:
        The maximum number of scatterings points to be saved for every photon hitting a DOM.
        Only the most recent positions are saved, older positions are overwritten if
        the maximum size is reached.
    :param UseI3PropagatorService:
        Use PROPOSAL and cmc to propagate initial particles (i.e. with NaN
        lengths) through the detector.
    :param UseGeant4:
        Enabling this setting will disable all cascade and muon light yield
        parameterizations. All particles will sent to Geant4 for a full
        simulation. This does **not** apply to muons that do have a length
        assigned. These are assumed to have been generated by MMC and
        their light is generated according to the usual parameterization.
    :param ParticleHistory:
        Store secondary particles produced by particle propagators (e.g.
        Geant4, PROPOSAL) in the MCTree.
    :param ParticleHistoryGranularity:
        When storing history in the MCTree, coalesce secondary EM cascades that
        are less than this distance apart.
    :param CrossoverEnergyEM:
        If set it defines the crossover energy between full Geant4 simulations and 
        light yield parameterizations for electro magnetic cascades. This only works
        when UseGeant4 is set to true. It works in conjunction with CrossoverEnergyHadron.
        If one of both is set to a positiv value greater 0 (GeV), the hybrid simulation
        is used.
        If CrossoverEnergyEM is set to None while CrossoverEnergyHadron is set so
        hybrid mode is working, GEANT4 is used for EM cascades.
        If CrossoverEnergyEM is set to 0 (GeV) while CrossoverEnergyHadron is set
        so hybrid mode is working, leptons and EM cascades will use parameterizations
        for the whole energy range.
    :param CrossoverEnergyHadron:
        If set it defines the crossover energy between full Geant4 simulations and
        light yield parameterizations for hadronic cascades. This only works when
        UseGeant4 is set to true. It works in conjunction with CrossoverEnergyEM.
        If one of both is set to a positiv value greater 0 (GeV), the hybrid simulation
        is used.
        If CrossoverEnergyHadron is set to None while CrossoverEnergyEM is set so
        hybrid mode is working, GEANT4 is used for hadronic cascades.
        If CrossoverEnergyHadron is set to 0 (GeV) while CrossoverEnergyHadron is
        set so hybrid mode is working, hadronic cascades will use parameterizations
        for the whole energy range.
    :param UseCascadeExtension:
    	If set, the cascade light emission parameterizations will include 
    	longitudinal extension. Otherwise, parameterized cascades will be 
    	treated as point-like. 
    :param DoNotParallelize:
        Try only using a single work item in parallel when running the
        OpenCL simulation. This might be useful if you want to run jobs
        in parallel on a batch system. This will only affect CPUs and
        will be a no-op for GPUs.
    :param DOMOversizeFactor:
        Set the DOM oversize factor. To disable oversizing, set this to 1.
    :param UnshadowedFraction:
        Fraction of photocathode available to receive light (e.g. unshadowed by the cable)
    :param HoleIceParameterization:
        Set this to a hole ice parameterization file. The default file contains the 
        coefficients for nominal angular acceptance correction due to hole ice (ice-models 
        project is required). Use file $I3_BUILD/ice-models/resources/models/angsens/as.nominal 
        for no hole ice parameterization.
    :param WavelengthAcceptance:
        If specified, use this wavelength acceptance to scale the generated
        Cherenkov spectrum rather than using the DOM acceptance modified for
        oversizing and angular acceptance.
    :param DOMRadius:
        Allow the DOMRadius to be set externally, for things like mDOMs.
    :param CableOrientation:
        Path to cable orientation file, e.g. $I3_BUILD/ice-models/resources/models/cable_position/orientation.led7.txt.
        If set, blocks photons that would have to pass through the best-fit
        position of the cable to reach a DOM. This reduces the isotropic
        efficiency by ~10%. Set to None to disable simulation of the cable
        shadow.
    :param OverrideApproximateNumberOfWorkItems:
        Allows to override the auto-detection for the maximum number of parallel work items.
        You should only change this if you know what you are doing.
    :param If:
        Python function to use as conditional execution test for segment modules.        
    :returns: the dictionary of keyword arguments passed to I3CLSimClientModule
    """

    from icecube import icetray, dataclasses, phys_services, clsim

    # warn the user in case they might have done something they probably don't want
    if UnWeightedPhotons and (DOMOversizeFactor != 1.):
        print("********************")
        print("Enabling the clsim.I3CLSimMakeHits() \"UnWeightedPhotons=True\" option without setting")
        print("\"DOMOversizeFactor=1.\" will still apply a constant weighting factor of DOMOversizeFactor**2.")
        print("If this is what you want, you can safely ignore this warning.")
        print("********************")
        
    if UnshadowedFraction<=0:
        raise RuntimeError("UnshadowedFraction must be a positive number")

    clsimParams = setupDetector(
        GCDFile=GCDFile,
        SimulateFlashers=bool(FlasherInfoVectName or FlasherPulseSeriesName),
        IceModelLocation=IceModelLocation,
        DisableTilt=DisableTilt,
        UnWeightedPhotons=UnWeightedPhotons,
        UnWeightedPhotonsScalingFactor=UnWeightedPhotonsScalingFactor,
        UseI3PropagatorService=UseI3PropagatorService,
        UseGeant4=UseGeant4,
        CrossoverEnergyEM=CrossoverEnergyEM,
        CrossoverEnergyHadron=CrossoverEnergyHadron,
        UseCascadeExtension=UseCascadeExtension,
        StopDetectedPhotons=StopDetectedPhotons,
        DOMOversizeFactor=DOMOversizeFactor,
        UnshadowedFraction=UnshadowedFraction,
        HoleIceParameterization=HoleIceParameterization,
        WavelengthAcceptance=WavelengthAcceptance,
        DOMRadius=DOMRadius,
        CableOrientation=CableOrientation,
        IgnoreSubdetectors=IgnoreSubdetectors,
    )

    converters = setupPropagators(RandomService, clsimParams,
        UseGPUs=UseGPUs,
        UseCPUs=UseCPUs,
        OverrideApproximateNumberOfWorkItems=OverrideApproximateNumberOfWorkItems,
        DoNotParallelize=DoNotParallelize,
        UseOnlyDeviceNumber=UseOnlyDeviceNumber,
        EnableDoubleBuffering=EnableDoubleBuffering,
        DoublePrecision=DoublePrecision,
        UseCUDA=UseCUDA,
    )
    # bind to a random port on localhost
    server = clsim.I3CLSimServer('tcp://127.0.0.1:*', clsim.I3CLSimStepToPhotonConverterSeries(converters))
    address = server.GetAddress()
    logging.log_info("Server listening at {}".format(address), unit="clsim")
    # stash server instance in the context to keep it alive
    tray.context[name+'CLSimServer'] = server

    if UseGPUs:
        if UseI3PropagatorService:
            logging.log_warn("Propagating muons and photons in the same process. This may starve your GPU.", unit="clsim")
        if UseGeant4:
            logging.log_warn("Running Geant and photon propagation in the same process. This will likely starve your GPU.", unit="clsim")

    module_config = \
    tray.Add(I3CLSimMakePhotonsWithServer, name,
        ServerAddress=address,
        DetectorSettings=clsimParams,
        MCTreeName=MCTreeName,
        OutputMCTreeName=OutputMCTreeName,
        FlasherInfoVectName=FlasherInfoVectName,
        FlasherPulseSeriesName=FlasherPulseSeriesName,
        PhotonSeriesName=PhotonSeriesName,
        MCPESeriesName=MCPESeriesName,
        RandomService=RandomService,
        ParticleHistory=ParticleHistory,
        ParticleHistoryGranularity=ParticleHistoryGranularity,
        ExtraArgumentsToI3CLSimClientModule=ExtraArgumentsToI3CLSimClientModule,
        If=If,
    )

    class GatherStatistics(icetray.I3Module):
        """Mimick the summary stage of I3CLSimModule::Finish()"""
        def Finish(self):
            if not 'I3SummaryService' in self.context:
                return
            summary = self.context['I3SummaryService']
            server = self.context[name+'CLSimServer']
            prefix = 'I3CLSimModule_'+name+'_makePhotons_clsim_'
            for k, v in server.GetStatistics().items():
                summary[prefix+k] = v
    tray.Add(GatherStatistics)

    return module_config

@my_traysegment
def I3CLSimMakePhotonsWithServer(tray, name,
                       ServerAddress,
                       DetectorSettings,
                       MCTreeName="I3MCTree",
                       OutputMCTreeName=None,
                       FlasherInfoVectName=None,
                       FlasherPulseSeriesName=None,
                       PhotonSeriesName="PhotonSeriesMap",
                       MCPESeriesName="MCPESeriesMap",
                       RandomService=None,
                       ParticleHistory=False,
                       ParticleHistoryGranularity=20*icetray.I3Units.m,
                       ExtraArgumentsToI3CLSimClientModule=dict(),
                       If=lambda f: True
                       ):
    """
    Propagate particles and photons up to PE level.
    Reads its particles from an I3MCTree and writes either an
    I3CompressedPhotonSeriesMap, an I3MCPESeriesMap, or both.

    Photon propagation is delegated to the I3CLSimServer listening at
    ServerAddress. This server may be shared by multiple client processes to
    maximize GPU utilization.

    If MMCTrackListName is set, this segment will assume that MMC has been
    applied to the I3MCTree and that MMC was *NOT* run using the "-recc"
    option.

    :param ServerAddress:
        The 0MQ address of an I3CLSimServer 
    :param DetectorSettings:
        The output of clsim.traysegments.common.setupDetector. These should be
        the same as those used to configure the server.
    :param MCTreeName:
        The name of the I3MCTree containing the particles to propagate.
    :param OutputMCTreeName:
        A copy of the (possibly sliced) MCTree will be stored as this name.
    :param FlasherInfoVectName:
        Set this to the name of I3FlasherInfoVect objects in the frame to
        enable flasher simulation. The module will read I3FlasherInfoVect objects
        and generate photons according to assumed parameterizations.
    :param FlasherPulseSeriesName:
        Set this to the name of an I3CLSimFlasherPulseSeries object in the frame to
        enable flasher/Standard Candle simulation.
        This cannot be used at the same time as FlasherInfoVectName.
        (I3CLSimFlasherPulseSeries objects are clsim's internal flasher
        representation, if "FlasherInfoVectName" is used, the I3FlasherInfoVect
        objects are converted to I3CLSimFlasherPulseSeries objects.)
    :param PhotonSeriesName:
        Configure this to enable writing an I3CompressedPhotonSeriesMap containing
        all photons that reached the DOM surface. If this is set to None,
        photons will be converted to PE immediately, drastically reducing
        memory overhead for bright events.
    :param MCPESeriesName:
        Configure this to enable writing an I3MCPESeriesMap.
    :param RandomService:
        Set this to an instance of a I3RandomService. Alternatively,
        you can specify the name of a configured I3RandomServiceFactory
        added to I3Tray using tray.AddService(). If you don't configure
        this, the default I3RandomServiceFactory will be used.
    :param ParticleHistory:
        Store secondary particles produced by particle propagators (e.g.
        Geant4, PROPOSAL) in the MCTree.
    :param ParticleHistoryGranularity:
        When storing history in the MCTree, coalesce secondary EM cascades that
        are less than this distance apart.
    :param If:
        Python function to use as conditional execution test for segment modules.        
    :returns: the dictionary of keyword arguments passed to I3CLSimClientModule
    """
    from icecube import icetray, dataclasses, phys_services, clsim

    # make sure the geometry is updated to the new granular format (in case it is supported)
    if hasattr(dataclasses, "I3ModuleGeo"):
        tray.AddModule("I3GeometryDecomposer", name + "_decomposeGeometry",
                       If=lambda frame: If(frame) and ("I3OMGeoMap" not in frame) and ("I3ModuleGeoMap" not in frame))

    if MCTreeName is None or MCTreeName=="":
        clSimMCTreeName=""
        if ChopMuons:
            raise RuntimeError("You cannot have \"MMCTrackListName\" enabled with no MCTree!")
    else:
        clSimMCTreeName=MCTreeName

    if FlasherInfoVectName is None or FlasherInfoVectName=="":
        if (FlasherPulseSeriesName is not None) and (FlasherPulseSeriesName!=""):
            SimulateFlashers=True
            clSimFlasherPulseSeriesName = FlasherPulseSeriesName
            clSimOMKeyMaskName = ""
        else:
            SimulateFlashers=False
            clSimFlasherPulseSeriesName = ""
            clSimOMKeyMaskName = ""
    else:
        if (FlasherPulseSeriesName is not None) and (FlasherPulseSeriesName!=""):
            raise RuntimeError("You cannot use the FlasherPulseSeriesName and FlasherInfoVectName parameters at the same time!")
        
        SimulateFlashers=True
        clSimFlasherPulseSeriesName = FlasherInfoVectName + "_pulses"
        clSimOMKeyMaskName = FlasherInfoVectName + "_OMKeys"
        
        tray.AddModule(clsim.FlasherInfoVectToFlasherPulseSeriesConverter,
                       name + "_FlasherInfoVectToFlasherPulseSeriesConverter",
                       FlasherInfoVectName = FlasherInfoVectName,
                       FlasherPulseSeriesName = clSimFlasherPulseSeriesName,
                       FlasherOMKeyVectName = clSimOMKeyMaskName,
                       If=If)

    if (OutputMCTreeName is not None) and (OutputMCTreeName != ""):
        # copy the MCTree to the requested output name
        def copyMCTree(frame, inputName, outputName, If_=None):
            if If_ is not None:
                if not If_(frame): return
            frame[outputName] = frame[inputName]
        tray.AddModule(copyMCTree, name + "_copyMCTree",
                       inputName=clSimMCTreeName,
                       outputName=OutputMCTreeName,
                       Streams=[icetray.I3Frame.DAQ],
                       If_=If)
        clSimMCTreeName = OutputMCTreeName
    else:
        clSimMCTreeName = clSimMCTreeName

    clsimModuleArgs = {
        'ServerAddress': ServerAddress,
        'PhotonSeriesMapName': PhotonSeriesName,
        'MCTreeName': clSimMCTreeName,
        'FlasherPulseSeriesName': clSimFlasherPulseSeriesName,
        'IgnoreSubdetectors': DetectorSettings['IgnoreSubdetectors'],
        'PhotonSeriesMapName': PhotonSeriesName or '',
        'MCPESeriesMapName': MCPESeriesName or '',
        'If': If,
    }
    clsimModuleArgs.update(**ExtraArgumentsToI3CLSimClientModule)

    # Set up light source -> step conversion if not already configured
    # This allows step generators with expensive initialization to be reused
    if not 'StepGenerator' in clsimModuleArgs:
        stepGenerator = clsim.I3CLSimLightSourceToStepConverterAsync(1)
        stepGenerator.SetLightSourceParameterizationSeries(DetectorSettings['ParameterizationList'])
        stepGenerator.SetMediumProperties(DetectorSettings['MediumProperties'])
        stepGenerator.SetRandomService(RandomService)
        stepGenerator.SetWlenBias(DetectorSettings['WavelengthGenerationBias'])
        propagators = []
        if DetectorSettings['UseI3PropagatorService']:
            from icecube.simprod.segments.PropagateMuons import make_standard_propagators
            pmap = make_standard_propagators(EmitTrackSegments=True, SplitSubPeVCascades=False)
            propagators.append(clsim.I3CLSimLightSourcePropagatorFromI3PropagatorService(pmap,ParticleHistory,ParticleHistoryGranularity))

        if DetectorSettings['UseGeant4']:
            propagators.append(clsim.I3CLSimLightSourcePropagatorGeant4(collectParticleHistory=ParticleHistory))

        stepGenerator.SetPropagators(propagators);
        clsimModuleArgs['StepGenerator'] = stepGenerator

    if clsimModuleArgs['MCPESeriesMapName']:
        # convert photons to MCPE in the worker thread of I3CLSimClientModule
        clsimModuleArgs['MCPEGenerator'] = clsim.I3CLSimPhotonToMCPEConverterForDOMs(
            RandomService,
            DetectorSettings['WavelengthAcceptance'],
            DetectorSettings['AngularAcceptance']
        )

    tray.AddModule('I3CLSimClientModule', name+"_makePhotons",
        **clsimModuleArgs
    )

    if MCPESeriesName:
        class GatherStatistics(icetray.I3Module):
            """Mimick the summary stage of I3PhotonToMCPEConverter::Finish()"""
            def __init__(self, ctx):
                icetray.I3Module.__init__(self, ctx)
                self.nhits = 0
            def DAQ(self, frame):
                self.nhits += sum((h.npe for hits in frame[MCPESeriesName].values() for h in hits))
                self.PushFrame(frame)
            def Finish(self):
                if not 'I3SummaryService' in self.context:
                    return
                summary = self.context['I3SummaryService']
                prefix = 'I3PhotonToMCPEConverter_'+name+'_makeHitsFromPhotons_clsim_make_hits_'
                summary[prefix+'NumGeneratedHits'] = self.nhits
        tray.Add(GatherStatistics)

    return clsimModuleArgs
