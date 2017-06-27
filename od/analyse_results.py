#!/home/nmrsu/anaconda3/bin/python
"""Read in the current experimental results, update the parameter estimates and
select the optimal point for the next observation.

The working directory should be specified as a command line argument. This
directory must contain the following input files:

    adapt_ccr_integrals: experimental observations of peak intensities
    adapt_ccr_tau: previously measured evolution times
    adapt_ccr_phase: previously measured phases (in degrees)
    adapt_ccr_theta: current best-estimate of parameter vector
    adapt_ccr_omega: resonance offsets (in angular units, s-1)

Following execution, two output files will be produced:

    adapt_ccr_newtau: next optimised sampling point
    adapt_ccr_newphase: next optimised sampling phase (in degrees)
    adapt_ccr_newtheta: updated estimate of parameter vector
    # adapt_ccr_dispersioncurve: dispersion curve for selection of new control point

Example (N.B. final slash is essential):
    >> ./analyse_results.py /opt/topspin3.5pl2/data/jc/nmr/chris_adaptive_methyl_ccr_080216/

"""

#from __future__ import division

import sys
import json
import numpy as np
from scipy import optimize, linalg
import pandas as pd

import response_surface_reduced as rs
import optimal_design_reduced as od
import information_matrix_reduced as im

# parse command line arguments
working_directory = sys.argv[1]
print("Working directory: " + working_directory)

def loadvar(name):
    # Load parameter from json file in working directory.
    f = open(working_directory + 'adapt_ccr_' + name, 'r')
    x = json.load(f)
    f.close()
    return x

def savevar(x, name):
    # Save parameter x to json file in working directory.
    f = open(working_directory + 'adapt_ccr_' + name, 'w')
    json.dump(x, f)
    f.close()

# load in integrals, tau, theta, omega from json files...
print "analyse_results.py - Input files:"
omega = np.array(loadvar('omega')).reshape((1,-1))
print omega
yobs = np.array(loadvar('integrals'))
yobs = yobs / yobs.max()
print yobs
taus = np.array(loadvar('taus'))#.reshape((-1,1))
print taus
#phases = np.array(loadvar('phases')).reshape((-1,1)) * np.pi / 180 # convert to rad
phases = np.array(loadvar('phases')) * np.pi / 180 # convert to rad
print phases
theta_0 = np.array(loadvar('theta0')).ravel()
print theta_0

# calculate current estimate of parameters
theta_hat, sigma, pcov = od.estimate_theta(yobs, taus, phases, omega, theta_0)
print "analyse_results.py - Current parameter estimate:"
#print np.column_stack((theta_hat, sigma))
relerr = 100 * sigma / theta_hat
print pd.DataFrame( \
            zip(theta_hat.tolist(), sigma.tolist(), relerr.tolist()), \
            columns=['fit','err','rel. err (%)'])

# get next point
(new_tau, new_phase) = od.direct_search_next_point(yobs, taus, phases, omega, theta_hat)
new_phase = int(np.round(new_phase * 180 / np.pi)) # convert to degrees

# save results
theta_hat = theta_hat.reshape((3,-1)).tolist()

savevar(new_tau, 'newtau')
savevar(new_phase, 'newphase')
savevar(theta_hat, 'newtheta')




