#!/usr/bin/env python


from icecube import icetray, dataclasses, simclasses, phys_services

rnd = phys_services.I3GSLRandomService(333)

def our_fn(x):
    print("Yes I was called with ", x, " but I don't know what that means.  Are you my daddy?")
    return x

gen = simclasses.I3SumGenerator(rnd, our_fn, 0, 10, 10, 5, 0.01, 100, 0.99, 100)

print("result:", gen.Generate(100))
