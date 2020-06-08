#!/usr/bin/env python

import unittest

from copy import deepcopy

from I3Tray import I3Tray, I3Units
from icecube import icetray
from icecube import dataclasses
from icecube import sim_services

class Source(icetray.I3Module):
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")

    def Configure(self):
        pass

    def Process(self):
        frame = icetray.I3Frame(icetray.I3Frame.DAQ)

        mctree = dataclasses.I3MCTree()

        ########################################################
        # This muon should be removed by I3InIceCORSIKATrimmer #
        ########################################################
        trimmed_muon = dataclasses.I3Particle()
        trimmed_muon.energy = 500 * I3Units.GeV
        trimmed_muon.type = dataclasses.I3Particle.MuMinus
        trimmed_muon.location_type = dataclasses.I3Particle.InIce

        
        ########################################################
        # This muon should be kept by I3InIceCORSIKATrimmer    #
        ########################################################
        kept_muon = dataclasses.I3Particle()
        kept_muon.energy = 1. * I3Units.TeV
        kept_muon.type = dataclasses.I3Particle.MuMinus
        kept_muon.location_type = dataclasses.I3Particle.InIce
        kept_muon.dir = dataclasses.I3Direction(0,0)

        ########################################################
        # This neutrino can be kept or dropped                 #
        ########################################################
        nu_mu = dataclasses.I3Particle()
        nu_mu.type = dataclasses.I3Particle.NuMu
        nu_mu.location_type = dataclasses.I3Particle.InIce

        ########################################################
        # This neutrino can be kept or dropped                 #
        ########################################################
        nu_e = dataclasses.I3Particle()
        nu_e.type = dataclasses.I3Particle.NuE
        nu_e.location_type = dataclasses.I3Particle.InIce

        ########################################################
        # This neutrino can be kept or dropped                 #
        ########################################################
        nu_tau = dataclasses.I3Particle()
        nu_tau.type = dataclasses.I3Particle.NuTau
        nu_tau.location_type = dataclasses.I3Particle.InIce

        ########################################################
        # I3InIceCORSIKATrimmer claims it's removing newly     #
        # childless parents, but really it's any childless     #
        # parent that's not InIce.                             #
        ########################################################
        mu_plus = dataclasses.I3Particle()
        mu_plus.type = dataclasses.I3Particle.MuPlus
        mu_plus.energy = 1. * I3Units.TeV
        mu_plus.location_type = dataclasses.I3Particle.InIce
        mu_plus.dir = dataclasses.I3Direction(0,0)

        mctree.add_primary(trimmed_muon)
        mctree.add_primary(kept_muon)
        mctree.add_primary(nu_mu)
        mctree.add_primary(nu_e)
        mctree.add_primary(nu_tau)
        mctree.append_child(nu_tau, mu_plus)

        print(mctree)
        
        frame["I3MCTree"] = mctree
        self.PushFrame(frame)


class TestInIceCORSIKATrimmer(unittest.TestCase):

    def setUp(self):
        self.tray = I3Tray()
        self.tray.AddModule(Source)

    def test_inice_corsika_trimmer_default(self):

        self.tray.AddModule("I3InIceCORSIKATrimmer")

        def TestModule(frame):
            self.assertTrue("I3MCTree" in frame)
            # started with 6 and now there should be 5
            # only the trimmed_muon below threshold should be trimmed
            print(len(frame["I3MCTree"]))
            mctree = frame["I3MCTree"]
            print(mctree)
            self.assertEqual(len(frame["I3MCTree"]), 5)

        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])
        self.tray.Execute(1)

    def test_inice_corsika_trimmer_drop_nu(self):

        self.tray.AddModule("I3InIceCORSIKATrimmer", RemoveNeutrinos = True)

        def TestModule(frame):
            self.assertTrue("I3MCTree" in frame)
            print(len(frame["I3MCTree"]))
            mctree = frame["I3MCTree"]
            print(mctree)
            # started with 6 and now there should be 1 
            self.assertEqual(len(frame["I3MCTree"]), 1)

        self.tray.AddModule(TestModule, streams = [icetray.I3Frame.DAQ])
        self.tray.Execute(1)

if __name__ == '__main__':
    unittest.main()
