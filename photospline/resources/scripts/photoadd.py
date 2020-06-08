#!/usr/bin/env python

"""
Add identical photonics tables together.

Usage: photoadd.py tab1.pt [tab2.pt tab3.pt ...] added.pt
"""

from icecube.photospline.photonics import Table, FITSTable
import sys

infiles, outfile = sys.argv[1:-1], sys.argv[-1]

if outfile.endswith('.fits'):
	TableClass = FITSTable
else:
	TableClass = Table

TableClass.stack(outfile, *infiles)

