#!/home/ucl/Enthought/Canopy_64bit/User/bin/python
"""Extract and analyse the final results of a methyl CCR measurement.

The working directory should be specified as a command line argument. This
directory must contain the following input files:

    adapt_ccr_tau
    adapt_ccr_integrals
    adapt_ccr_theta
    adapt_ccr_omega

Example (N.B. final slash is essential):
    >> ./final_results.py /opt/topspin3.2/data/jc/nmr/chris_adaptive_methyl_ccr_090116/

"""

from __future__ import division

import sys
import json
import numpy as np
from math import pi
from scipy import optimize
import matplotlib.pyplot as plt
import pandas as pd

import response_surface

# parse command line arguments
working_directory = sys.argv[1]
print("Working directory: " + working_directory)

def loadvar(name):
    # Load parameter from json file in working directory.
    f = open(working_directory + 'adapt_ccr_' + name, 'r')
    x = json.load(f)
    f.close()
    return x

tau = np.array(loadvar('tau')).reshape((-1,1))
integrals = np.array(loadvar('integrals'))
theta = np.array(loadvar('theta'))
omega = np.array(loadvar('omega'))

x = pd.DataFrame(integrals,index=tau)
with pd.option_context('display.max_rows',999):
    print "\nINTEGRALS:\n"
    print x
    print

plt.hist(tau,bins=20)
plt.show()

# normalise integrals by their maximum value
yobs = integrals / integrals.max()

def resid(theta):
    return (response_surface.y(tau,theta,omega) - yobs).ravel()

out = optimize.leastsq(resid, theta, full_output=1)
theta_hat = out[0]
# calculate uncertainties
# https://mail.scipy.org/pipermail/scipy-user/2013-March/034316.html
pcov = out[1]

chi2 = sum(resid(theta_hat)**2) / (yobs.size - theta_hat.size)

print("chi2 = %g\n\n" % chi2)

sigma = np.sqrt(pcov.diagonal() * chi2)

relerr=100. * sigma / theta_hat
fit_results = pd.DataFrame(zip(theta_hat.tolist(), sigma.tolist(), relerr.tolist()), columns=['fit','err','rel. err (%)'])
print("FIT RESULTS:\n")
print fit_results
