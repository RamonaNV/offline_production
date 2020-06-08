#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example script for MuonGun

IceCube is illuminated isotropically with muons. This is not useful for
production but can be used to calculate effective areas.

"""
import icecube
import icecube.icetray
import icecube.dataclasses
import icecube.dataio
import icecube.phys_services
import icecube.MuonGun
import icecube.MuonGun.segments
from I3Tray import I3Tray


def main(output_file="./floodlight.i3.gz", num_events=100,
         gamma_index=2, energy_offset=1e3, energy_min=1e4, energy_max=1e7,
         coszen_min=-1., coszen_max=1.):
    tray = I3Tray()

    try:
        tray.context["I3RandomService"] =\
            icecube.phys_services.I3SPRNGRandomService(1, 10000, 1)
    except AttributeError:
        tray.context["I3RandomService"] =\
            icecube.phys_services.I3GSLRandomService(1)

    surface = icecube.MuonGun.Cylinder(
        length=1600.*icecube.icetray.I3Units.m,
        radius=800.*icecube.icetray.I3Units.m)

    spectrum = icecube.MuonGun.OffsetPowerLaw(
        gamma_index,
        energy_offset*icecube.icetray.I3Units.GeV,
        energy_min*icecube.icetray.I3Units.GeV,
        energy_max*icecube.icetray.I3Units.GeV)

    # Illuminate IceCube isotropically with muons.
    generator = icecube.MuonGun.Floodlight(
        surface, spectrum, coszen_min, coszen_max)

    tray.AddSegment(icecube.MuonGun.segments.GenerateBundles, "muongun",
                    Generator=generator,
                    NEvents=num_events,
                    GCDFile="")

    # Calculate effective area.
    tray.AddModule(effective_area, "effective_area",
                   generator=generator,
                   Streams=[icecube.icetray.I3Frame.DAQ])

    tray.AddModule("I3Writer", "writer",
                   filename=output_file,
                   Streams=[icecube.icetray.I3Frame.DAQ])

    

    tray.Execute()
    

    return


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
    return True


if __name__ == "__main__":
    import optparse

    import numpy

    parser = optparse.OptionParser("%prog [options]\n\n" + __doc__)

    parser.add_option("-o", "--out",
                      dest="output_file",
                      default="floodlight.i3.gz",
                      help="path to output i3-file ['%default']")
    parser.add_option("-n", "--num",
                      dest="num_events",
                      default=100,
                      type="int",
                      help="number of events [%default]")
    parser.add_option("--gamma",
                      dest="gamma_index",
                      default=2.,
                      type="float",
                      help="spectral index [%default]")
    parser.add_option("--offset",
                      dest="energy_offset",
                      default=1e3,
                      type="float",
                      help="offset energy in GeV [%default]")
    parser.add_option("--emin",
                      dest="energy_min",
                      default=1e4,
                      type="float",
                      help="minimum energy in GeV [%default]")
    parser.add_option("--emax",
                      dest="energy_max",
                      default=1e7,
                      type="float",
                      help="maximum energy in GeV [%default]")
    parser.add_option("--tmin",
                      dest="zenith_min",
                      default=0.,
                      type="float",
                      help="minimum zenith angle in deg [%default]")
    parser.add_option("--tmax",
                      dest="zenith_max",
                      default=180.,
                      type="float",
                      help="maximum zenith angle in deg [%default]")

    opts, args = parser.parse_args()

    main(opts.output_file, opts.num_events, opts.gamma_index,
         opts.energy_offset, opts.energy_min, opts.energy_max,
         numpy.cos(numpy.deg2rad(opts.zenith_max)),
         numpy.cos(numpy.deg2rad(opts.zenith_min)))
