"""Compute information matrices."""

import numpy as np
from numba import jit

import response_surface as rs


@jit(nopython=True)
def Mdiscrete(t, phi, theta, omega):
    """
    Compute information matrix for a discrete experiment with control points
    (t, phi) conditioned on the parameter vector, theta, and the offsets, omega.

    M =  SUM  F(t,phi; theta,omega) . F(t,phi; theta,omega)'
        expts

    where F is the Jacobian matrix.

    Args: (where M = number of time points, N = number of peaks)
        t: (M x 1) array of evolution times
        phi: (M x 1) array of phases (in radians)
        theta: flattened (3 x N) array (with array.ravel).
                theta[0,:] = amplitudes
                theta[1,:] = lambda (auto-relaxation rate)
                theta[2,:] = S2tc (in ns)
        omega: (1 x N) array of resonance offsets, in angular units (s-1).

    Returns:
        y: (M x N) array of intensities
    """

    #if phi is None:
    #    phi = np.zeros_like(t)

#    M = np.zeros( (theta.size,theta.size) );
    M = np.zeros( (len(theta), len(theta)) );

    for i in range(len(t)):
        F = rs.jac(t[i], phi[i], theta, omega)
        M = M + np.dot(F,F.T)
#        M = M + F.dot(F.T)
#    for tau, phase in np.nditer([t, phi]):
#        F = rs.jac(tau,phase,theta,omega)
#        M = M + F.dot(F.T)

    return M
