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
# gamma = 40. # resolution enhancement
# alpha = 2 # 1 = max. sensitivty; 2 = constant noise

# File names:
schedule_file = 'vclist' # location of vclist file containing NUWS schedule

################################################################################

# little function to find closest entry in list
def closest(array, X):
    idx = [(np.abs(array - x)).argmin() for x in X]
    return idx

# Calculate required apodization function:
t = t0 + np.linspace(0, in_f*npoints, npoints) # evolution times
if nuws_type == 'Jdec':
    s = np.cos(np.pi * Jcc * t) # expected signal envelope
    ns = (1. / s)**2
    # calculate zero crossings (up to 4th = 100 ms for 35 Hz)
    nulls = np.array([0.5,1.5,2.5,3.5]) / Jcc
    nulls = nulls[nulls < max(t)]
    nullidx = closest(t, nulls)
    deadpoints = []
    for idx in nullidx:
        deadpoints.extend(range(idx-dx,idx+dx))
    ns[deadpoints] = 1 # set points near zero crossings to nulls
elif nuws_type == 'resolution':
    ns = np.exp(alpha * gamma * t)
ns = np.rint(ns) # round ns to integer

# Write schedule file to disk
with open(schedule_file,'w') as f:
    for n in ns:
        f.write("{:d}\n".format(int(n)))

print('NUWS schedule successfully written to {:s}.'.format(schedule_file))
print('Mean sampling density: {}.\n'.format(np.mean(ns)))

