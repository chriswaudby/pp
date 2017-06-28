#!/home/nmrsu/anaconda3/bin/python
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

import sys
import json
import numpy as np
from math import pi
from scipy import optimize
import matplotlib.pyplot as plt

import response_surface as rs
import optimal_design as od
import utils

N_seed = 5

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


# parse command line arguments
working_directory = sys.argv[1]
print("Working directory: " + working_directory)


# load data
yobs = np.array(loadvar('integrals'))
yobs = yobs / yobs.max() * 10.
omega = np.array(loadvar('omega'))
taus = np.array(loadvar('taus'))#.reshape((-1,1))
phases = np.array(loadvar('phases')) * np.pi / 180 # convert to rad
theta_0 = np.array(loadvar('theta0')).ravel()
N_expts = len(taus)
N_spins = len(omega)

#utils.pprint_array(omega, label='omega')
#utils.pprint_array(yobs, label='integrals')
#utils.pprint_array(taus, label='taus')
#utils.pprint_array(phases, label='phases')
#utils.pprint_array(theta_0, 'theta0')


# calculate running estimates of parameters
S2tc = []
S2tc_err = []
R2 = []
R2_err = []
theta_hat = theta_0
for i in range(N_expts):
    if i<N_seed:
        R2.append(theta_0[N_spins:2*N_spins])
        R2_err.append( 0. * R2[0] )
        S2tc.append(theta_0[2*N_spins:3*N_spins])
        S2tc_err.append( 0. * S2tc[0] )
    else:
        #theta_hat, sigma = od.estimate_theta(yobs[:i,:], taus[:i], phases[:i], omega, theta_0)
        theta_hat, sigma = od.estimate_theta(yobs[:i,:], taus[:i], phases[:i], omega, theta_hat)
        R2.append(theta_hat[N_spins:2*N_spins])
        R2_err.append( sigma[N_spins:2*N_spins] )
        S2tc.append(theta_hat[2*N_spins:3*N_spins])
        S2tc_err.append( sigma[2*N_spins:3*N_spins] )
R2 = np.array(R2)
R2_err = np.array(R2_err)
S2tc = np.array(S2tc)
S2tc_err = np.array(S2tc_err)

print("Running parameter estimates:")
utils.pprint_array(R2, label='R2(running)')
utils.pprint_array(S2tc, label='S2tc(running)')

print("Final parameter estimates:")
print_estimate(theta_hat, sigma)

for i in range(N_spins):
    plt.subplot(2,3,i+1)
    #plt.errorbar(np.arange(N_expts), S2tc[:,i], yerr=S2tc_err[:,i], fmt='ob-')
    plt.fill_between(np.arange(N_expts), R2[:,i]-R2_err[:,i], R2[:,i]+R2_err[:,i], facecolor='palegreen')
    plt.plot(np.arange(N_expts), R2[:,i], 'g-')
    plt.fill_between(np.arange(N_expts), S2tc[:,i]-S2tc_err[:,i], S2tc[:,i]+S2tc_err[:,i], facecolor='mistyrose')
    plt.plot(np.arange(N_expts), S2tc[:,i], 'r-')
    plt.xlabel('Iteration number')
    #plt.ylabel('S2tc estimate / ns')
    plt.ylim(ymin=0)
    #plt.title('spin {:g}'.format(i+1))
plt.show()

#for i in range(N_spins):
#    plt.subplot(2,3,i+1)
#    plt.fill_between(np.arange(N_expts), R2[:,i]-R2_err[:,i], R2[:,i]+R2_err[:,i], facecolor='honeydew')
#    plt.plot(np.arange(N_expts), R2[:,i], 'g-')
#    plt.xlabel('Iteration number')
#    plt.ylabel('R2(13C) estimate / s-1')
#    plt.ylim(ymin=0)
#    #plt.title('spin {:g}'.format(i+1))
#plt.show()

# plot distribution of sampled times
#plt.hist(taus,bins=40)
#plt.show()



