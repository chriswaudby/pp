#!/usr/bin/env python
import numpy as np

# Acqusition parameters:
npoints = 128 # number of complex points (td / 2)
in_f= 315.6e-6 # increment for indirect evolution period (in seconds)
t0 = in_f / 2. # evolution period for first point (e.g. half-dwell)
# Note: for J deconvolution, t0 should include any additional periods during
# which the chemical shift is refocused but the homonuclear coupling continues
# to evolve.

# Window function parameters: select Jdec/resolution as required
# 1. J decoupling:
nuws_type = 'Jdec'
Jcc = 35. # homonuclear coupling constant for decoupling
dx = 7 # number of complex points to omit either side of nulls
# 2. Resolution enhancement:
# nuws_type = 'resolution'
# gamma = 10. # resolution enhancement

# File names:
schedule_file = 'vclist' # location of vclist file containing NUWS schedule
input_file = 'raw.ft1' # nmrPipe format file, before processing of indirect dimension
output_file = 'scaled.ft1' # output file name
nuslist_file = 'nuslist' # nuslist output file name (for homonuclear J decoupling)

################################################################################

# little function to find closest entry in list
def closest(array, X):
    idx = [(np.abs(array - x)).argmin() for x in X]
    return idx

# Load input data and NUWS schedule
with open(input_file, 'rb') as fid:
    header = np.fromfile(fid, np.int8, count=2048)
    y = np.fromfile(fid, np.float32)
y=y.reshape((2*npoints,-1))
ns = np.loadtxt(schedule_file)

# Calculate required apodization function:
t = t0 + np.linspace(0, in_f*npoints, npoints) # evolution times
if nuws_type == 'Jdec':
    s = np.cos(np.pi * Jcc * t) # expected signal envelope
    window = 1. / s 
elif nuws_type == 'resolution':
    window = np.exp(gamma * t)
h = window / ns # window = ns * h

# Form scaled output
yscaled = 0. * y
yscaled[0::2,:] = y[0::2,:] * h[:, None] # re
yscaled[1::2,:] = y[1::2,:] * h[:, None] # im

# Write output
yscaled = yscaled.reshape((-1))
with open(output_file, 'wb') as fid:
    header.tofile(fid)
    yscaled.tofile(fid)

# Generate nuslist for J decoupling:
if nuws_type == 'Jdec':
    # calculate zero crossings (up to 4th = 100 ms for 35 Hz)
    nulls = np.array([0.5,1.5,2.5,3.5]) / Jcc
    nulls = nulls[nulls < max(t)]
    nullidx = closest(t, nulls)
    deadpoints = []
    for idx in nullidx:
        deadpoints.extend(range(idx-dx,idx+dx))
    with open(nuslist_file,'w') as f:
        for i in range(npoints):
            if i not in deadpoints:
                f.write("{:d}\n".format(i))
