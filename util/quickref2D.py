#!/usr/bin/env python

from acqus import get_entry

import argparse

# parse command line arguments
parser = argparse.ArgumentParser(description='Extract a parameter from bruker acqus format file')
parser.add_argument('nuc',choices=['N','C'],help='heteronucleus')
parser.add_argument('pipe_shift',type=float,help='nmrPipe H2O shift estimate')
args = parser.parse_args()
xpipe = args.pipe_shift
nuc = args.nuc

dx = xpipe - 4.7
o1 = float(get_entry('acqus','O1')) 
bf1 = float(get_entry('acqus','BF1'))
xcar = o1/bf1 + dx
sf1 = float(get_entry('acqus','SFO1'))
v0 = sf1 / (1+xcar*1e-6)
sf2 = float(get_entry('acqu2s','SFO1'))
if nuc=='C':
    xi = 0.25144953
else:
    xi = 0.101329118
v2 = xi*v0
ycar = (sf2-v2)/v2*1e6

print('Quick referencing:')
print('XCAR = {:.3f} ppm'.format(xcar))
print('YCAR = {:.3f} ppm'.format(ycar))
