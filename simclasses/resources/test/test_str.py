#!/usr/bin/env python
import unittest

def make_instance(className):
    import icecube.simclasses
    instance = eval("icecube.simclasses." + className + "()")
    return instance


def is_standard_repr(s):
    import re
    return re.search("object at [01]x[0-9a-f]+", s) is not None


class TestStr(unittest.TestCase):

    stringList = []  # to store produced strings

    def test_empty(self):
        # skipping maps which don't have proper __str__ method
        # skipping I3MMCTrack*, where __str__ is not implemented yet
        classNameList = ['CorsikaLongProfile',
                         'CorsikaLongStep',
                         'I3CorsikaShowerInfo',
                         'I3MCPE',
                         'I3MCPESeries',
                         'I3MCPulse',
                         'I3MCPulseSeries']

        nBad = 0
        for className in classNameList:
            instance = make_instance(className)
            s = str(instance)
            self.stringList.append("<" + className+"(empty)>\n"+s)
            if is_standard_repr(s):
                nBad += 1

        self.assertEqual(nBad, 0)

    def test_CorsikaLongProfile(self):
        from icecube.simclasses import CorsikaLongStep, CorsikaLongProfile

        s = CorsikaLongStep()
        s.depth = 10
        s.numCharged = 20
        p = CorsikaLongProfile()
        p.append(s)
        p.append(s)

        s = str(p)
        self.stringList.append("<CorsikaLongProfile>\n"+s)
        self.assertFalse(is_standard_repr(s))

    def test_CorsikaShowerInfo(self):
        from icecube.simclasses import CorsikaLongStep, CorsikaLongProfile, I3CorsikaShowerInfo

        i = I3CorsikaShowerInfo()
        s = CorsikaLongStep()
        s.depth = 10
        s.numCharged = 20
        p = i.longProfile
        p.append(s)
        p.append(s)

        s = str(i)
        self.stringList.append("<I3CorsikaShowerInfo>\n"+s)
        self.assertFalse(is_standard_repr(s))

    def test_I3MCPESeries(self):
        from icecube.simclasses import I3MCPE, I3MCPESeries

        x = I3MCPE()
        seq = I3MCPESeries()
        seq.append(x)
        seq.append(x)

        s = str(seq)
        self.stringList.append("<I3MCPESeries>\n"+s)
        self.assertFalse(is_standard_repr(s))

    def test_I3MCPulseSeries(self):
        from icecube.simclasses import I3MCPulse, I3MCPulseSeries

        x = I3MCPulse()
        seq = I3MCPulseSeries()
        seq.append(x)
        seq.append(x)

        s = str(seq)
        self.stringList.append("<I3MCPulseSeries>\n"+s)
        self.assertFalse(is_standard_repr(s))

    @classmethod
    def tearDownClass(self):
        # print the produced strings, test outputs are in arbitary order
        from sys import stderr
        stderr.flush()
        stderr.write("\nTestOutput\n----------\n" + "\n".join(self.stringList) + "\n")
        stderr.flush()

if __name__ == "__main__":
    unittest.main()
