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

import sys
import json
import numpy as np
#from scipy import optimize, linalg
#import pandas as pd

import response_surface as rs
import optimal_design as od
import information_matrix as im
import utils
import parameters

max_aq_time = 0.2

# parse command line arguments
working_directory = sys.argv[1]
print("Working directory: " + working_directory)


def print_estimate(theta_hat, sigma, label='ESTIMATE'):
    n_spins = len(theta_hat) // 3
    for spin in range(n_spins):
        print('{:s}\t{:d}\tamp\t{:.3f}\t{:.3f}'.format(label, spin+1,theta_hat[spin],sigma[spin]))
    for spin in range(n_spins):
        print('{:s}\t{:d}\tR2\t{:.3f}\t{:.3f}'.format(label, spin+1,theta_hat[spin+n_spins],sigma[spin+n_spins]))
    for spin in range(n_spins):
        print('{:s}\t{:d}\tS2tc\t{:.3f}\t{:.3f}'.format(label, spin+1,theta_hat[spin+2*n_spins],sigma[spin+2*n_spins]))



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
print("analyse_results.py - Input files:")

omega = np.array(loadvar('omega')).reshape((1,-1))

yobs = np.array(loadvar('integrals'))
yobs = yobs / yobs.max() * 10.
taus = np.array(loadvar('taus'))#.reshape((-1,1))
phases = np.array(loadvar('phases')) * np.pi / 180 # convert to rad
theta_0 = np.array(loadvar('theta0')).ravel()

# select just first spin for debugging
#yobs = yobs[:,0]
#omega = omega[:,0:2]
#theta_0 = theta_0[::5]

#utils.pprint_array(omega, label='omega')
#utils.pprint_array(yobs, label='integrals')
#utils.pprint_array(taus, label='taus')
#utils.pprint_array(phases, label='phases')
#utils.pprint_array(theta_0, 'theta0')


# calculate current estimate of parameters
theta_hat, sigma = od.estimate_theta(yobs, taus, phases, omega, theta_0)
print("analyse_results.py - Current parameter estimate:")
print_estimate(theta_hat, sigma)

# get next point
(new_tau, new_phase) = od.direct_search_next_point(yobs, taus, phases, omega, theta_hat, tau_max=max_aq_time)
new_phase = int(np.round(new_phase * 180 / np.pi)) # convert to degrees

print("analyse_results.py - Optimal acquistion time = {:.6f} s, phase = {:.1f}".format(new_tau, new_phase))

# save results
theta_hat = theta_hat.reshape((3,-1)).tolist()

savevar(new_tau, 'newtau')
savevar(new_phase, 'newphase')
savevar(theta_hat, 'newtheta')




