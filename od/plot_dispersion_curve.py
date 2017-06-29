#!/home/ucl/Enthought/Canopy_64bit/User/bin/python
"""Plot the previously calculated dispersion curve for adaptive sampling of methyl CCR.

The working directory should be specified as a command line argument. This
directory must contain the following input files:

    adapt_ccr_dispersioncurve

Example:
    >> ./plot_dispersion_curve.py /opt/topspin3.2/data/jc/nmr/chris_150116

"""

from __future__ import division

import sys
import json
import numpy as np
from math import pi
import matplotlib.pyplot as plt

# parse command line arguments
working_directory = sys.argv[1]
print("Working directory: " + working_directory)

def loadvar(name):
    # Load parameter from json file in working directory.
    f = open(working_directory + 'adapt_ccr_' + name, 'r')
    x = json.load(f)
    f.close()
    return x

dispersioncurve = loadvar('dispersioncurve')
t, m = zip(*dispersioncurve)

plt.plot(t,m)
plt.show()
