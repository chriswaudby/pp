#!/usr/bin/python

import sys
import argparse
import numpy as np
import nmrglue as ng

parser = argparse.ArgumentParser(description='Sum a collection of bruker ser files. Data will be truncated to fit shortest ser file.')
parser.add_argument('-out', metavar='outputser', dest='out', help='Output filename', required=True)
parser.add_argument('serfiles', metavar='inputser', type=argparse.FileType('rb'), nargs='+',
                    help='ser files to be summed')
args = parser.parse_args()

fileformat = '<i4' # little-endian

# load first file
serfile = args.serfiles[0]
print('Loading {0}'.format(serfile))
data = np.frombuffer(serfile.read(), dtype=fileformat)

if len(args.serfiles)>1:
    for serfile in args.serfiles[1:]:
        print('Adding {0}'.format(serfile))
        data1 = np.frombuffer(serfile.read(), dtype=fileformat)
        maxlength = min([data.size, data1.size])
        if data.size != data1.size:
            print("*** WARNING! ser files are being truncated! ***")
        data = data[:maxlength] + data1[:maxlength]

# write output
print('Writing output to {0}'.format(args.out))
with open(args.out,'wb') as f:
    f.write(data.astype(fileformat).tostring())

print('Finished!')
