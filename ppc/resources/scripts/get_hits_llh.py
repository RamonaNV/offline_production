#!/usr/bin/env python
"""
Extract hits from pulse series to text file for llh
"""
import argparse
import os
import itertools
import numpy as np
from icecube import icetray, dataio, dataclasses, phys_services, TopologicalSplitter
from I3Tray import I3Tray, I3Units

curr = 0

def extract(frame, pulses, residual, window, outdir, subdir, ini):
    global curr

    def write_tws(fme, key, fobj):
        if fme.Has(key):
            fmeobj = fme[key]
            for dom in fmeobj.keys():
                for tw in fmeobj[dom]:
                    f.write('{} {} {} {}\n'.format(dom.string, dom.om, tw.start, tw.stop))
        
    all_pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(
        frame, pulses)

    if subdir is None:
        outdir = os.path.join(outdir,
                              '{}.{}'.format(frame['I3EventHeader'].run_id,
                                            frame['I3EventHeader'].event_id))
    else:
        outdir = os.path.join(outdir, subdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    stop = np.inf
    if residual is not None:
        med = np.median([p.time for p in itertools.chain(*all_pulses.values())])
        stop = med+residual

    with open(os.path.join(outdir, 'dat{:03d}'.format(curr)), 'w') as f:
        for dom in all_pulses.iterkeys():
            pulses = all_pulses[dom]
            for pulse in pulses:
                if stop-window < pulse.time < stop:                    
                    f.write('{} {} {} {}\n'.format(dom.string, dom.om, pulse.time, pulse.charge))

    with open(os.path.join(outdir, 'ert{:03d}'.format(curr)), 'w') as f:
        write_tws(frame, 'CalibrationErrata', f)
        write_tws(frame, 'SaturationWindows', f)        

    if ini is not None and frame.Has(ini):
        p = frame[ini]
        energy = 1e5 if np.isnan(p.energy) else p.energy
        with open(os.path.join(outdir, 'ini{:03d}'.format(curr)), 'w') as f:
            f.write(' '.join(map(str,
                                 [p.pos.x, p.pos.y, p.pos.z,
                                 p.dir.zenith/I3Units.deg, p.dir.azimuth/I3Units.deg,
                                 energy, p.time]))+'\n')
            if p.is_track:
                f.write('\n')
                f.write('0 10 3 0')

    curr += 1


def main():
    parser = argparse.ArgumentParser(
        description='This program extracts hits from a pulse series into a text file')
    parser.add_argument('inp', nargs='+', help='input files')
    parser.add_argument('-g', '--gcd', type=str, default=None, help='GCD file')
    parser.add_argument('-o', '--outdir', type=str, default='dats',
                        help='output directory')
    parser.add_argument('-P', '--preserve', default=False, action='store_true',
                        help='use the input filename as part of the outdir instead of eventid.runid')
    parser.add_argument('-p', '--pulses', type=str, default=None,
                        help='specify the pulse series to process, defaults to output of TopologicalSplitter')
    parser.add_argument('-r', '--residual', type=float, default=None,
                        help='specify the residual time after the median of pulse-series after which pulses are not saved. Default is to save all hits.')
    parser.add_argument('-w', '--window', type=float, default=np.inf,
                        help='specify the time window length before the med+residual stop time over which pulses are saved. Default is infinite.')
    parser.add_argument('--raw', type=str,
                        default='InIcePulses', help='specify the Q-frame pulse to process with TopologicalSplitter')
    parser.add_argument('--ini', type=str,
                        default=None, help='specify key in frame to use as seed for extraction to ini')

    args = parser.parse_args()

    flist = sorted(args.inp)
    for fname in flist:
        tray = I3Tray()
        if args.gcd is None:
            tray.Add('I3Reader', Filename=fname)
        else:
            tray.Add('I3Reader', Filenamelist=[args.gcd, fname])

        pulses = args.pulses
        if pulses is None:
            pulses = 'SplitTopoPulses'
            tray.Add('I3TopologicalSplitter', InputName=args.raw,
                     OutputName=pulses,
                     Multiplicity=5,
                     TimeWindow=2000*I3Units.ns,
                     TimeCone=500*I3Units.ns,
                     XYDist=300*I3Units.m,
                 ZDomDist=15)

        if args.preserve:
            subdir = os.path.splitext(os.path.basename(fname))[0]
        else:
            subdir = None

        tray.Add(extract,
                 pulses=pulses,
                 residual=args.residual,
                 window=args.window,
                 outdir=args.outdir,
                 subdir=subdir,
                 ini=args.ini,
                 If=lambda frame:frame.Has(pulses))

        tray.Execute()


if __name__ == '__main__':
    main()
