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

tau = np.array([1.1, 3, 6, 10, 15, 25, 40])

lines = lines[18:]
n = 0
for line in lines:
    n += 1
    plt.subplot(4,6,n)
    dat = line.strip().split()
    pid = dat[0]
    I = float(dat[17])
    err = float(dat[18])

    ratio = np.array([])
    ratio_err = np.array([])

    for i in range(len(tau)):
        y1 = float(dat[7+25+i]) * I / 24
        y1_err = float(dat[7+25+i]) * err / 24
        y2 = float(dat[25+i]) * I / 16
        y2_err = float(dat[25+i]) * err / 16
        ratio = np.append(ratio, y1/y2)
        ratio_err = np.append(ratio_err, y1/y2 * np.sqrt((y1_err/y1)**2 + (y2_err/y2)**2))
        #plt.errorbar(tau[i], y1, yerr=y1_err, fmt='ob')
        #plt.errorbar(tau[i], y2, yerr=y2_err, fmt='sr')
    plt.errorbar(tau, ratio, yerr=ratio_err, fmt='o-')
    plt.ylim((0,1))
    if n==24:
        plt.show()
        n = 0
plt.show()

