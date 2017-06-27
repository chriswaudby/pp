"""Define a response surface for methyl CCR.
The CCR rates sigma and eta, and relaxation during INEPT transfers, are all calculated from S2tc."""

#from __future__ import division   # needed for python 2

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
wC = 2*np.pi * gC/gH * 800e6   # 700 MHz
P2cosb = -1./3.  # for 109 degree tetrahedral geometry
carbon_csa = 18e-6   # 18 ppm for Ile, 25 ppm for Leu/Val

# pre-factors for calculation relaxation rates
c_CHCH = 2/45. * (mu0*hbar*gC*gH / (4*pi*rCH**3))**2 * 1e-9
c_CHC = 2/5. * (mu0/(4*pi)) * P2cosb * rCH**-3 * hbar * gH * gC * wC * carbon_csa * 1e-9
c_HHHC = 1/5. * (mu0/(4*pi))**2 * rCH**-3 * rHH**-3 * hbar**2 * gH**3 * gC * 1e-9
c_HHHH = 9/20. * (mu0*hbar*gH**2 / (4*pi*rHH**3))**2 * 1e-9

# test functions
def calc_sigma(S2tc):
    return c_CHCH * S2tc

def calc_eta(S2tc):
    return c_CHC * S2tc

# methyl 1J(CH) scalar coupling, in Hz
J = 125.
piJ = pi * J
TAU = 0.5 / J

@jit
def calc_DELTA(S2tc):
    """Compute relaxation during INEPT transfers due to CH/CH and CH/HH dipole/dipole CCR (in macromolecular limit)."""
    return np.exp(-TAU * c_HHHH * S2tc) * np.cosh(TAU * c_HHHC * S2tc);

@jit
def y(tau,phi,theta,omega):
    """Compute expected signal intensities.

    y = A * exp(-lambda t) * (
         Iout exp([-3 sigma - 2 eta]*t) cos([omega + 3 pi J]*t + phi)
       + Iin exp([sigma - 2/3 eta]*t) cos([omega + pi J]*t + phi)
       + Iin exp([sigma + 2/3 eta]*t) cos([omega - pi J]*t + phi)
       + Iout exp([-3 sigma + 2 eta]*t) cos([omega - 3 pi J]*t + phi)
         )
    
    Iin = 3 + 3*DELTA(S2tc)
    Iout = 3 - DELTA(S2tc)
    
    sigma = c_CHCH * tau_c
    eta = c_CHC * tau_c

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
    ph = phi.reshape((-1,1))

    parameter_matrix = theta.reshape((3,-1))
    A = parameter_matrix[0,:].reshape((1,-1))
    lam = parameter_matrix[1,:].reshape((1,-1))
    S2tc = parameter_matrix[2,:].reshape((1,-1))

#    print('theta:')
#    print(theta)
#    print('parameter_matrix')
#    print(parameter_matrix)
#    print('amplitudes:')
#    print(A)
#    print('lambda:')
#    print(lam)
#    print('S2tc:')
#    print(S2tc)
#    print('omega:')
#    print(omega)
#    print('tau:')
#    print(tau)
#    print('t:')
#    print(t)

    sigma = c_CHCH * S2tc
    eta = c_CHC * S2tc

    DELTA = calc_DELTA(S2tc)
    Iouter = 3 + 3*DELTA
    Iinner = 3 - DELTA

#    print('sigma')
#    print(sigma)
#    print('eta')
#    print(eta)
#    print('c_HHHC*S2tc')
#    print(c_HHHC*S2tc)

#    print('DELTA')
#    print(DELTA)
#    print('Iinner')
#    print(Iinner)
#    print('Iouter')
#    print(Iouter)

#    print( np.exp(-lam*t) )
#    print(Iouter*np.exp((-3*sigma - 2*eta)*t)*np.cos((omega + 3*piJ)*t + ph))
#    print(Iinner*np.exp((sigma - 2*eta/3)*t)*np.cos((omega + piJ)*t + ph))

    return A * np.exp(-lam*t) * ( \
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
    for n in range(N):
        A = parameter_matrix[0,n]
        lam = parameter_matrix[1,n]
        S2tc = parameter_matrix[2,n] 

        sigma = c_CHCH * S2tc
        eta = c_CHC * S2tc
        DELTA = calc_DELTA(S2tc)
        Iinner = 3 + 3*DELTA
        Iouter = 3 - DELTA

        w = omega[0,n]

        exp_lam = np.exp(-lam*t)
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
               (TAU*c_CHCH*DELTA + (DELTA-3)*(3*c_CHCH+2*c_CHC)*t - TAU*c_HHHC*np.exp(-TAU*c_CHCH*S2tc)*np.sinh(TAU*c_HHHC*S2tc) ) * \
               np.cos((w + 3*piJ)*t + phi) \
            + np.exp((sigma - 2*eta/3)*t) * \
               (-3*TAU*c_CHCH*DELTA + (DELTA+1)*(3*c_CHCH-2*c_CHC)*t + 3*TAU*c_HHHC*np.exp(-TAU*c_CHCH*S2tc)*np.sinh(TAU*c_HHHC*S2tc) ) * \
               np.cos((w + piJ)*t + phi) \
            + np.exp((sigma + 2*eta/3)*t) * \
               (-3*TAU*c_CHCH*DELTA + (DELTA+1)*(3*c_CHCH+2*c_CHC)*t + 3*TAU*c_HHHC*np.exp(-TAU*c_CHCH*S2tc)*np.sinh(TAU*c_HHHC*S2tc) ) * \
               np.cos((w - piJ)*t + phi) \
            + np.exp((-3*sigma + 2*eta)*t) * \
               (TAU*c_CHCH*DELTA + (DELTA-3)*(3*c_CHCH-2*c_CHC)*t - TAU*c_HHHC*np.exp(-TAU*c_CHCH*S2tc)*np.sinh(TAU*c_HHHC*S2tc) ) * \
               np.cos((w - 3*piJ)*t + phi) )

    return F
