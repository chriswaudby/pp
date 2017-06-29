"""Define a reduced response surface for methyl CCR.
The CCR rates sigma, eta and deltaR2 are all calculated from S2tc."""

from __future__ import division

import numpy as np
from numba import jit
from math import pi

# physical constants
mu0 = 4e-7 * np.pi
hbar = 1.055e-34
gH = 2.675e8
gC = 6.726e7
rCH = 1.117e-10 # Tugarinov 2004
rHH = np.sqrt(3) * rCH * np.sin(110.4*np.pi/180) # Tugarinov 2004
wC = 2*np.pi * 176e6
P2cosb = -.228e30 * rCH**3
carbon_csa = 18e-6

c_sigma = (2/45)*(mu0*hbar*gH*gC/(4*np.pi*rCH**3))**2  # CH/CH dipole/dipole
c_eta = -(2/5)*(mu0*hbar*gH*gC/(4*np.pi*rCH**3)) * wC * carbon_csa * P2cosb  # CH/C dipole/CSA
c_zeta = (9/20)*(mu0*hbar*gH**2/(4*np.pi*rHH**3))**2  # HH/HH dipole/dipole

# methyl 1J(CH) scalar coupling, in Hz
J = 125.
piJ = pi * J
DELTA = 0.5 / J

@jit
def calc_sigma(S2tc):
    """Compute CH/CH dipole/dipole CCR rate (in macromolecular limit)."""
    return c_sigma * S2tc

@jit
def calc_eta(S2tc):
    """Compute CH/C dipole/CSA CCR rate (in macromolecular limit)."""
    return c_eta * S2tc

@jit
def calc_zeta(S2tc):
    """Compute HH/HH dipole/dipole CCR rate (in macromolecular limit)."""
    return c_zeta * S2tc

@jit
def y(tau,phi,theta,omega):
    """Compute expected signal intensities.

    y = (A/128) * exp(-lambda t) * (
            Iout exp([-3 sigma - 2 eta]*t) cos([omega + 3 pi J]*t + phi)
            + Iin exp([sigma - 2/3 eta]*t) cos([omega + pi J]*t + phi)
            + Iin exp([sigma + 2/3 eta]*t) cos([omega - pi J]*t + phi)
            + Iout exp([-3 sigma + 2 eta]*t) cos([omega - 3 pi J]*t + phi)
            )
    
    Iin = 9R**2 - 4R + 11
    Iout = 9R**2 + 24R + 15
    R = exp(-DELTA * zeta)
    
    sigma = c_sigma * tau_c
    eta = c_eta * tau_c
    zeta = c_zeta * tau_c

    Args: (where M = number of time points, N = number of peaks)
        t: (M x 1) array of evolution times
        phi: (M x 1) array of phases (in radians)
        theta: flattened (3 x N) array (with array.ravel)
                theta[0,:] = amplitudes
                theta[1,:] = lambda (auto-relaxation rate)
                theta[2,:] = S2tc (in ns)
        omega: (1 x N) array of resonance offsets, in angular units (s-1)

    Returns:
        y: (M x N) array of intensities
    """

    t = tau.reshape((-1,1))
    ph = -1 * phi.reshape((-1,1))

    parameter_matrix = theta.reshape((3,-1))
    A = parameter_matrix[0,:].reshape((1,-1))
    lam = parameter_matrix[1,:].reshape((1,-1))
    S2tc = parameter_matrix[2,:].reshape((1,-1))
    
    sigma = calc_sigma(S2tc*1e-9)
    eta = calc_eta(S2tc*1e-9)
    zeta = calc_zeta(S2tc*1e-9)

    R = np.exp(-DELTA*zeta)

    Iinner = 9*R**2 - 4*R + 11
    Iouter = 9*R**2 + 24*R + 15

    return (A/128) * np.exp(-lam*t) * ( \
        Iouter*np.exp((-3*sigma - 2*eta)*t)*np.cos((omega + 3*piJ)*t + ph) \
        + Iinner*np.exp((sigma - 2*eta/3)*t)*np.cos((omega + piJ)*t + ph) \
        + Iinner*np.exp((sigma + 2*eta/3)*t)*np.cos((omega - piJ)*t + ph) \
        + Iouter*np.exp((-3*sigma + 2*eta)*t)*np.cos((omega - 3*piJ)*t + ph))


@jit(nopython=True)
def jac(t,phi,theta,omega):
    """Compute Jacobian of response surface (at a given evolution time).

    F_ij = dy_j / dq_i

    Args: (where N = number of peaks)
        t (float): evolution time
        phi (float): phase (in radians)
        theta: flattened (3 x N) array (with array.ravel)
                theta[0,:] = amplitudes
                theta[1,:] = lambda (auto-relaxation rate)
                theta[2,:] = S2tc (in ns)
        omega: (1 x N) array of resonance offsets, in angular units (s-1)

    Returns:
        F: (3N x N) array with jacobian matrix
    """

    N = omega.size
    F = np.zeros((3*N,N), dtype=np.float64)

    parameter_matrix = theta.reshape((3,-1))
    # loop over spins:
    for n in xrange(N):
        A = parameter_matrix[0,n]
        lam = parameter_matrix[1,n]
        S2tc = parameter_matrix[2,n] * 1e-9
        
        sigma = c_sigma * S2tc
        eta = c_eta * S2tc
        zeta = c_zeta * S2tc

        R = np.exp(-DELTA*zeta)
        Iinner = 9*R**2 - 4*R + 11
        Iouter = 9*R**2 + 24*R + 15
        
        w = omega[0,n]

        exp_lam = np.exp(-lam*t) / 128
        x1 = np.exp((-3*sigma - 2*eta)*t) * np.cos((w + 3*piJ)*t + phi)
        x2 = np.exp((sigma - 2*eta/3)*t) * np.cos((w + piJ)*t + phi)
        x3 = np.exp((sigma + 2*eta/3)*t) * np.cos((w - piJ)*t + phi)
        x4 = np.exp((-3*sigma + 2*eta)*t) * np.cos((w - 3*piJ)*t + phi)

        # dy/dA:
        F[3*n,n] = exp_lam * ( \
            Iouter*x1 \
            + Iinner*x2 \
            + Iinner*x3 \
            + Iouter*x4)
        # dy/dlambda:
        F[3*n+1,n] = -t * A * exp_lam * ( \
            Iouter*x1 \
            + Iinner*x2 \
            + Iinner*x3 \
            + Iouter*x4)
        # dy/d(S2tc)
        F[3*n+2,n] = A * exp_lam * ( \
            np.exp((-3*sigma - 2*eta)*t) * \
               (Iouter*(-3*c_sigma*t - 2*c_eta*t) - 24*c_zeta*R*DELTA - 18*c_zeta*R**2*DELTA) * \
               np.cos((w + 3*piJ)*t + phi) \
            + np.exp((sigma - 2*eta/3)*t) * \
               (Iinner*(c_sigma*t - 2*c_eta*t/3) + 4*c_zeta*R*DELTA - 18*c_zeta*R**2*DELTA) * \
               np.cos((w + piJ)*t + phi) \
            + np.exp((sigma + 2*eta/3)*t) * \
               (Iinner*(c_sigma*t + 2*c_eta*t/3) + 4*c_zeta*R*DELTA - 18*c_zeta*R**2*DELTA) * \
               np.cos((w - piJ)*t + phi) \
            + np.exp((-3*sigma + 2*eta)*t) * \
               (Iouter*(-3*c_sigma*t + 2*c_eta*t) - 24*c_zeta*R*DELTA - 18*c_zeta*R**2*DELTA) * \
               np.cos((w - 3*piJ)*t + phi) )

    return F
