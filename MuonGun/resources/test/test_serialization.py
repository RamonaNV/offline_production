#!/usr/bin/env python

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--read", default=None)
parser.add_argument("--write", default=None)
args = parser.parse_args()

from icecube import icetray, dataio, MuonGun, phys_services

generator = MuonGun.StaticSurfaceInjector()

def read(fname='muongun_serialization_test.i3'):
    f = dataio.I3File(fname)
    frame = f.pop_frame()
    f.close()

    newgenerator = frame['Generator']
    assert newgenerator.surface == generator.surface

def write(fname='muongun_serialization_test.i3'):

    frame = icetray.I3Frame()
    frame['Generator'] = generator
    
    f = dataio.I3File(fname, 'w')
    f.push(frame)
    f.close()

if args.read:
    read(args.read)
elif args.write:
    write(args.write)
else:
    write()
    read()



