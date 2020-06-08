"""
Tray segments for MuonGun simulations
"""
import math

import icecube.icetray
import icecube.dataclasses
import icecube.MuonGun
import icecube.MuonGun.segments


@icecube.icetray.traysegment
def GenerateCosmicRayMuons(tray, name,
                           mctree_name="I3MCTree_preMuonProp",
                           num_events=1,
                           flux_model="Hoerandel5_atmod12_SIBYLL",
                           gamma_index=2.,
                           energy_offset=700.,
                           energy_min=1e4,
                           energy_max=1e7,
                           cylinder_length=1600.,
                           cylinder_radius=800.,
                           cylinder_x=0.,
                           cylinder_y=0.,
                           cylinder_z=0.,
                           inner_cylinder_length=500.,
                           inner_cylinder_radius=150.,
                           inner_cylinder_x=46.3,
                           inner_cylinder_y=-34.9,
                           inner_cylinder_z=-300.,
                           use_inner_cylinder=False):
    r"""Generate atmospheric muons with MuonGun

    This segment generates atm. muons with MuonGun; it is intended to be
    used within :class:`icecube.simprod.modules.MuonGunGenerator`.

    .. note::
        Only single muons are generated, no bundles. For a effective
        generation, the in-ice muon scpectrum is approximated with
        :math:`dP/dE_{\mu} \approx (E_{\mu} + b)^{-\gamma}`.

    :param str mctree_name:
        Output name for :cpp:class:`I3MCTree`
    :param int num_events:
        Number of events to generate
    :param str flux_model:
        Name of primary cosmic-ray flux parametrization
    :param float gamma_index:
        Spectral index :math:`\gamma` of power-law approximation
    :param float energy_offset:
        Energy offset :math:`b` of power-law approximation
    :param float energy_min:
        Minimum generated energy in GeV
    :param float energy_max:
        Maximum generated energy in GeV
    :param float cylinder_length:
        Length of injection cylinder in m
    :param float cylinder_radius:
        Radius of injection cylinder in m
    :param float cylinder_x:
        Cartesian coordinates :math:`(x, y, z)` of injection cylinder's
        center in m
    :param float cylinder_y:
        Cartesian coordinates :math:`(x, y, z)` of injection cylinder's
        center in m
    :param float cylinder_z:
        Cartesian coordinates :math:`(x, y, z)` of injection cylinder's
        center in m
    :param float inner_cylinder_length:
        Length of inner injection cylinder in m; relevant for DeepCore
        simulations.
    :param float inner_cylinder_radius:
        Radius of inner injection cylinder in m; relevant for DeepCore
        simulations.
    :param float inner_cylinder_x:
        Cartesian coordinates :math:`(x, y, z)` of inner injection
        cylinder's center in m; relevant for DeepCore simulations.
    :param float inner_cylinder_y:
        Cartesian coordinates :math:`(x, y, z)` of inner injection
        cylinder's center in m; relevant for DeepCore simulations.
    :param float inner_cylinder_z:
        Cartesian coordinates :math:`(x, y, z)` of inner injection
        cylinder's center in m; relevant for DeepCore simulations.
    :param bool use_inner_cylinder:
        Use inner injection cylinder together with constant surface
        scaling function; relevant for DeepCore simulations.

    """
    # Default: use Hoerandel as a template for generating muons.
    model = icecube.MuonGun.load_model(flux_model)

    # Generate only single muons, no bundles.
    model.flux.max_multiplicity = 1

    # Default: cylinder aligned with z-axis at detector center as injection
    # surface.
    outsurface_center = icecube.dataclasses.I3Position(
        cylinder_x*icecube.icetray.I3Units.m,
        cylinder_y*icecube.icetray.I3Units.m,
        cylinder_z*icecube.icetray.I3Units.m)

    outsurface = icecube.MuonGun.Cylinder(
        length=cylinder_length*icecube.icetray.I3Units.m,
        radius=cylinder_radius*icecube.icetray.I3Units.m,
        center=outsurface_center)

    # Draw energies from a power law with offset.
    spectrum = icecube.MuonGun.OffsetPowerLaw(
        gamma=gamma_index,
        offset=energy_offset*icecube.icetray.I3Units.GeV,
        min=energy_min*icecube.icetray.I3Units.GeV,
        max=energy_max*icecube.icetray.I3Units.GeV)

    # Set up the generator. This gets stored in a special frame.
    if use_inner_cylinder:
        # Use energy-dependent scaling if generating atm. muon in DeepCore.
        insurface_center = icecube.dataclasses.I3Position(
            inner_cylinder_x*icecube.icetray.I3Units.m,
            inner_cylinder_y*icecube.icetray.I3Units.m,
            inner_cylinder_z*icecube.icetray.I3Units.m)

        insurface = icecube.MuonGun.Cylinder(
            length=inner_cylinder_length*icecube.icetray.I3Units.m,
            radius=inner_cylinder_radius*icecube.icetray.I3Units.m,
            center=insurface_center)

        scaling = icecube.MuonGun.ConstantSurfaceScalingFunction(insurface)

        generator = icecube.MuonGun.EnergyDependentSurfaceInjector(
            surface=outsurface, flux=model.flux, energy=spectrum,
            radius=model.radius, scaling=scaling)
    else:
        generator = icecube.MuonGun.StaticSurfaceInjector(
            outsurface, model.flux, spectrum, model.radius)

    tray.AddModule('I3MuonGun::GeneratorModule',name,
                   Generator=num_events*generator)

    tray.AddModule("I3MuonGun::WeightCalculatorModule", "%s_weights" % name,
                   Generator=generator,
                   Model=model)

    tray.AddModule("Rename", "%s_prepropMCTree" % name,
                   keys=["I3MCTree", mctree_name])

    return


@icecube.icetray.traysegment
def GenerateSingleMuons(tray, name,
                        Surface=None,
                        GCDFile=None,
                        GeometryMargin=60.*icecube.icetray.I3Units.m,
                        NumEvents=100,
                        FromEnergy=10.*icecube.icetray.I3Units.TeV,
                        ToEnergy=10.*icecube.icetray.I3Units.PeV,
                        BreakEnergy=1.*icecube.icetray.I3Units.TeV,
                        GammaIndex=2.,
                        ZenithRange=[0., 180.*icecube.icetray.I3Units.deg]):
    r"""Generate single muons with MuonGun

    This segment injects single muons with MuonGun isotropically. This
    useful for calculating effective areas.

    .. note::
        For a effective generation, the in-ice muon scpectrum is
        approximated with
        :math:`dP/dE_{\mu} \approx (E_{\mu} + b)^{-\gamma}`.

    :param icecube.MuonGun.SamplingSurface Surface:
        Injection surface; if not specified and if `GCDFile` is
        ``None``, a cylinder with a length of 1600 m and a radius of
        800 m is used.
    :param str GCDFile:
        Path to geometry file; if specified, the injection surface will
        be an :class:`icecube.MuonGun.ExtrudedPolygon` around the
        in-ice DOMs. Remote file locations are supported via
        :class:`icecube.dataio.I3FileStager`.
    :param float GeometryMargin:
        If `GCDFile` is specified, this margin (in meters) is used to
        expand the convex hull around the in-ice DOMs.
    :param int NumEvents:
        Number of events to generate
    :param float FromEnergy:
        Minimum generated energy in GeV
    :param float ToEnergy:
        Maximum generated energy in GeV
    :param float GammaIndex:
        Spectral index :math:`\gamma` of power-law approximation
    :param float BreakEnergy:
        Energy offset :math:`b` of power-law approximation
    :param list ZenithRange:
        Zenith angle range in radians for isotropic injection

    """
    if Surface is None and GCDFile is None:
        surface = icecube.MuonGun.Cylinder(
            length=1600.*icecube.icetray.I3Units.m,
            radius=800.*icecube.icetray.I3Units.m)
    elif Surface is not None and GCDFile is None:
        surface = Surface
    elif Surface is None and GCDFile is not None:
        if "I3FileStager" in tray.context:
            handle = tray.context["I3FileStager"].GetReadablePath(GCDFile)
            GCDFile = str(handle)
        surface = icecube.MuonGun.ExtrudedPolygon.from_file(
            GCDFile, padding=GeometryMargin*icecube.icetray.I3Units.m)
    else:
        icecube.icetray.logging.log_fatal(
            "Surface and GCDFile are mutually exclusive.",
            unit="GenerateSingleMuons")

    spectrum = icecube.MuonGun.OffsetPowerLaw(
        GammaIndex,
        BreakEnergy*icecube.icetray.I3Units.GeV,
        FromEnergy*icecube.icetray.I3Units.GeV,
        ToEnergy*icecube.icetray.I3Units.GeV)

    # Illuminate IceCube isotropically with muons.
    generator = NumEvents*icecube.MuonGun.Floodlight(
        surface, spectrum, math.cos(ZenithRange[1]), math.cos(ZenithRange[0]))

    tray.AddModule("I3MuonGun::GeneratorModule", name,
                   Generator=generator)

    # Calculate effective area.
    def effective_area(frame, generator):
        mctree = frame["I3MCTree"]
        primary = mctree.primaries[0]
        muon = mctree.get_daughters(primary)[0]
        bundle = icecube.MuonGun.BundleConfiguration(
            [icecube.MuonGun.BundleEntry(0, muon.energy)])
        fluence = generator.generated_events(primary, bundle)
        if fluence > 0.:
            area = 1./fluence
            frame["MCMuon"] = muon
            frame["MuonEffectiveArea"] = icecube.dataclasses.I3Double(area)
        else:
            icecube.icetray.logging.log_warn(
                "Fluence value of {0:f} encountered.".format(fluence),
                unit="GenerateSingleMuons")
        return True

    tray.Add(effective_area,
             generator=generator,
             Streams=[icecube.icetray.I3Frame.DAQ])

    tray.Add("Rename",
             keys=["I3MCTree", "I3MCTree_preMuonProp"])

    return

@icecube.icetray.traysegment
def GenerateNaturalRateMuons(tray, name,
                        mctree_name="I3MCTree_preMuonProp",
                        flux_model="GaisserH4a_atmod12_SIBYLL",
                        Surface=None,
                        GCDFile=None,
                        GeometryMargin=60.*icecube.icetray.I3Units.m,
                        NumEvents=100):
    r"""Generate single muons with MuonGun

    This segment injects single muons with MuonGun isotropically. This
    useful for calculating effective areas.

    .. note::
        For a effective generation, the in-ice muon scpectrum is
        approximated with
        :math:`dP/dE_{\mu} \approx (E_{\mu} + b)^{-\gamma}`.

    :param str mctree_name:
        Output name for :cpp:class:`I3MCTree`
    :param str flux_model:
        Name of primary cosmic-ray flux parametrization
    :param icecube.MuonGun.SamplingSurface Surface:
        Injection surface; if not specified and if `GCDFile` is
        ``None``, a cylinder with a length of 1600 m and a radius of
        800 m is used.
    :param str GCDFile:
        Path to geometry file; if specified, the injection surface will
        be an :class:`icecube.MuonGun.ExtrudedPolygon` around the
        in-ice DOMs. Remote file locations are supported via
        :class:`icecube.dataio.I3FileStager`.
    :param float GeometryMargin:
        If `GCDFile` is specified, this margin (in meters) is used to
        expand the convex hull around the in-ice DOMs.
    :param int NumEvents:
        Number of events to generate

    """
    if Surface is None and GCDFile is None:
        surface = icecube.MuonGun.Cylinder(
            length=1600.*icecube.icetray.I3Units.m,
            radius=800.*icecube.icetray.I3Units.m)
    elif Surface is not None and GCDFile is None:
        surface = Surface
    elif Surface is None and GCDFile is not None:
        if "I3FileStager" in tray.context:
            handle = tray.context["I3FileStager"].GetReadablePath(GCDFile)
            GCDFile = str(handle)
        surface = icecube.MuonGun.ExtrudedPolygon.from_file(
            GCDFile, padding=GeometryMargin*icecube.icetray.I3Units.m)
    else:
        icecube.icetray.logging.log_fatal(
            "Surface and GCDFile are mutually exclusive.",
            unit="GenerateSingleMuons")

    # Sample from given spectral flux model
    model = icecube.MuonGun.load_model(flux_model)

    flux = model.flux
    radii = model.radius
    energies = model.energy

    generator = NumEvents*icecube.MuonGun.NaturalRateInjector(surface, flux, energies)

    tray.AddModule("I3MuonGun::GeneratorModule", name,
                   Generator=generator)
    tray.AddModule('I3MuonGun::WeightCalculatorModule', 'weight',
                   Model=icecube.MuonGun.BundleModel(flux, radii, energies),
                   Generator=generator)


    # Calculate effective area.
    def effective_area(frame, generator):
        mctree = frame["I3MCTree"]
        primary = mctree.primaries[0]
        muon = mctree.get_daughters(primary)[0]
        bundle = icecube.MuonGun.BundleConfiguration(
            [icecube.MuonGun.BundleEntry(0, muon.energy)])
        fluence = generator.generated_events(primary, bundle)
        if fluence > 0.:
            area = 1./fluence
            frame["MCMuon"] = muon
            frame["MuonEffectiveArea"] = icecube.dataclasses.I3Double(area)
        else:
            icecube.icetray.logging.log_warn(
                "Fluence value of {0:f} encountered.".format(fluence),
                unit="GenerateSingleMuons")
        return True

    tray.Add(effective_area,
             generator=generator,
             Streams=[icecube.icetray.I3Frame.DAQ])

    tray.AddModule("Rename", "%s_prepropMCTree" % name,
                   keys=["I3MCTree", mctree_name])

    return
