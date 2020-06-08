#!/usr/bin/env python
import unittest

class TestEquality(unittest.TestCase):

    def test_CorsikaLongStep(self):
        from icecube.simclasses import CorsikaLongStep
        a = CorsikaLongStep()
        a.depth = 1
        a.numCharged = 2
        b = CorsikaLongStep()
        b.depth = 1
        b.numCharged = 2
        c = CorsikaLongStep()
        c.depth = 1
        c.numCharged = 4

        self.assertEqual(a, a)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)

    def test_CorsikaLongProfile(self):
        from icecube.simclasses import CorsikaLongProfile, CorsikaLongStep
        x = CorsikaLongStep()
        x.depth = 1
        x.numCharged = 2

        a = CorsikaLongProfile()
        a.append(x)
        b = CorsikaLongProfile()
        b.append(x)
        c = CorsikaLongProfile()
        c.append(x)
        c.append(x)

        self.assertEqual(a, a)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)

    def test_I3CorsikaShowerInfo(self):
        from icecube.simclasses import I3CorsikaShowerInfo, CorsikaLongProfile, CorsikaLongStep

        def make():
            # fill fields with some sensible numbers, avoid NaN
            x = I3CorsikaShowerInfo()
            x.firstIntHeight = 1.1
            x.firstIntDepth = 2.2
            x.obsLevelHeight = 3.3
            x.ghMaxNum = 4.4
            x.ghStartDepth = 5.5
            x.ghRedChiSqr = 6.6
            x.resampleRadius = 7.7
            x.ghMaxDepth = 8.8
            x.ghLambdaa = 9.9
            x.ghLambdab = 10.0
            x.ghLambdac = 11.1
            return x

        x = CorsikaLongStep()
        x.depth = 1
        x.numCharged = 2

        a = make()
        a.longProfile.append(x)
        b = make()
        b.longProfile.append(x)
        c = make()
        d = make()
        d.longProfile.append(x)
        d.longProfile.append(x)

        self.assertEqual(a, a)
        self.assertEqual(a, b)
        self.assertEqual(c, c)
        self.assertNotEqual(a, c)
        self.assertNotEqual(a, d)
        self.assertNotEqual(c, d)

    def test_I3MCPE(self):
        from icecube.simclasses import I3MCPE, I3MCPESeries

        def make():
            x = I3MCPE()
            x.time = 10.1
            x.npe = 3
            return x

        a = make()
        b = make()
        c = make()
        c.npe = 4

        self.assertEqual(a, a)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)

        d = I3MCPESeries()
        d.append(a)
        e = I3MCPESeries()
        e.append(b)
        f = I3MCPESeries()
        f.append(c)

        self.assertEqual(d, d)
        self.assertEqual(d, e)
        self.assertNotEqual(d, f)

    def test_I3MCPulse(self):
        from icecube.simclasses import I3MCPulse, I3MCPulseSeries

        def make():
            x = I3MCPulse()
            x.time = 10.1
            x.charge = 3.3
            x.source = I3MCPulse.PE
            return x

        a = make()
        b = make()
        c = make()
        c.charge = 4.4

        self.assertEqual(a, a)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)

        d = I3MCPulseSeries()
        d.append(a)
        e = I3MCPulseSeries()
        e.append(b)
        f = I3MCPulseSeries()
        f.append(c)

        self.assertEqual(d, d)
        self.assertEqual(d, e)
        self.assertNotEqual(d, f)

if __name__ == "__main__":
    unittest.main()
