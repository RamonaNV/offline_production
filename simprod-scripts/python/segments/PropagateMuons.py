"""Tray segments for muon propagation

"""
import os

import icecube
import icecube.icetray
import icecube.dataclasses
import icecube.phys_services
import icecube.sim_services
import icecube.simclasses
import icecube.cmc
import icecube.PROPOSAL
import json

default_media_definition = os.path.expandvars(
    "$I3_BUILD/PROPOSAL/resources/config_icesim.json")


@icecube.icetray.traysegment
def PropagateMuons(tray, name,
                   RandomService=None,
                   CylinderRadius=None,
                   CylinderLength=None,
                   SaveState=True,
                   InputMCTreeName="I3MCTree_preMuonProp",
                   OutputMCTreeName="I3MCTree",
                   **kwargs):
    r"""Propagate muons.

    This segment propagates muons through ice with ``PROPOSAL``; it
    simulates lepton decays and energy losses due to ionization,
    bremsstrahlung, photonuclear interactions, and pair production.

    :param I3RandomService RandomService:
        Random number generator service
    :param float CylinderRadius:
        Radius of the target volume in m
        (this param is now depricated, use the config file in the detector configuration)
    :param float CylinderLength:
        Full height of the target volume in m
        (this param is now depricated, use the config file in the detector configuration)
    :param bool SaveState:
        If set to `True`, store the state of the supplied RNG.
    :param str InputMCTree:
        Name of input :cpp:class:`I3MCTree` frame object
    :param str OutputMCTree:
        Name of output :cpp:class:`I3MCTree` frame object
    :param \**kwargs:
        Additional keyword arguments are passed to
        :func:`icecube.simprod.segments.make_propagator`.

    """
    if CylinderRadius is not None:
        icecube.icetray.logging.log_warn(
            "The CylinderRadius now should be set in the configuration file in the detector configuration")
    if CylinderLength is not None:
        icecube.icetray.logging.log_warn(
            "The CylinderLength now should be set in the configuration file in the detector configuration")
    propagators = make_standard_propagators(**kwargs)

    # Set up propagators.
    if "I3ParticleTypePropagatorServiceMap" in tray.context:
        propagator_map = tray.context["I3ParticleTypePropagatorServiceMap"]
        for k, v in propagators.items():
            propagator_map[k] = v
    else:
        propagator_map = propagators

    if SaveState:
        rng_state = InputMCTreeName+"_RNGState"
    else:
        rng_state = ""

    tray.AddModule("I3PropagatorModule", name+"_propagator",
                   PropagatorServices=propagator_map,
                   RandomService=RandomService,
                   InputMCTreeName=InputMCTreeName,
                   OutputMCTreeName=OutputMCTreeName,
                   RNGStateName=rng_state)

    # Add empty MMCTrackList objects for events that have none.
    def add_empty_tracklist(frame):
        if "MMCTrackList" not in frame:
            frame["MMCTrackList"] = icecube.simclasses.I3MMCTrackList()
        return True

    tray.AddModule(add_empty_tracklist, name+"_add_empty_tracklist",
                   Streams=[icecube.icetray.I3Frame.DAQ])

    return

def make_standard_propagators(SplitSubPeVCascades=True,
                              EmitTrackSegments=True,
                              MaxMuons=10,
                              PROPOSAL_config_file=default_media_definition):
    """
    Set up standard propagators (PROPOSAL for muons and taus, CMC for cascades)
    :param bool SplitSubPeVCascades:
        Split cascades into segments above 1 TeV. Otherwise, split only above 1 PeV.
    :param bool EmitTrackSegments:
        Emit constant-energy track slices in addition to stochastic losses
        (similar to the output of I3MuonSlicer)
    :param str PROPOSAL_config_file:
        Path to PROPOSAL config file
    Keyword arguments will be passed to I3PropagatorServicePROPOSAL
    """
    from icecube.icetray import I3Units

    cascade_propagator = icecube.cmc.I3CascadeMCService(
        icecube.phys_services.I3GSLRandomService(1))  # Dummy RNG
    cascade_propagator.SetEnergyThresholdSimulation(1*I3Units.PeV)
    if SplitSubPeVCascades:
        cascade_propagator.SetThresholdSplit(1*I3Units.TeV)
    else:
        cascade_propagator.SetThresholdSplit(1*I3Units.PeV)
    cascade_propagator.SetMaxMuons(MaxMuons)
    muon_propagator = icecube.PROPOSAL.I3PropagatorServicePROPOSAL(
            config_file=PROPOSAL_config_file, slice_tracks=EmitTrackSegments)
    propagator_map =\
        icecube.sim_services.I3ParticleTypePropagatorServiceMap()

    for pt in "MuMinus", "MuPlus", "TauMinus", "TauPlus":
        key = getattr(icecube.dataclasses.I3Particle.ParticleType, pt)
        propagator_map[key] = muon_propagator

    for pt in "DeltaE", "Brems", "PairProd", "NuclInt", "Hadrons",\
              "EMinus", "EPlus":
        key = getattr(icecube.dataclasses.I3Particle.ParticleType, pt)
        propagator_map[key] = cascade_propagator
    
    return propagator_map

