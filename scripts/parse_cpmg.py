#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

input_file = 'nlin.tab'

f = open(input_file, 'r')
lines = f.readlines()
f.close()

# 0 = id
# 17 = intensity
# 18 = error
# 25+ = data

ncyc = np.array([40, 1, 36, 2, 32, 3, 28, 4, 24, 5, 20, 6, 18, 7, 16, 8, 14, 9, 12, 10])
vCPMG = ncyc * 25

lines = lines[18:]
n = 0
for line in lines:
    n += 1
    plt.subplot(6,12,n)
    dat = line.strip().split()
    pid = dat[0]
    I = float(dat[17])
    err = float(dat[18])
    R2 = np.array([])
    R2err = np.array([])
    yref = float(dat[25]) * I
    yref_err = float(dat[25]) * err
    for i in range(len(ncyc)):
        y = float(dat[26+i]) * I
        y_err = float(dat[26+i]) * err
        R2 = np.append(R2,-1/.04 * np.log(y/yref))
        R2err = np.append(R2err,1/.04 * np.sqrt((yref_err/yref)**2 + (y_err/y)**2) )
    plt.errorbar(vCPMG, R2, yerr=R2err, fmt='o')
    plt.ylim((0,40))
    if n==12*6:
        plt.show()
        n = 0
plt.show()

