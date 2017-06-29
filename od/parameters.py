import numpy as np

import response_surface as rs
import utils



VERBOSE = False

n_spins = 3

noise = 0.1 * 16  # 16 = rough scaling based on (2*Iinner + 2*Iouter)

# set up 'true' parameters (these are based on real expt 119 measurement)
#omega = np.array([[8.89,  415.,  -674.]])
#amplitude = np.array([1.36, 1.46, 1.67])
#R2 = np.array([61.4, 62., 48.8])
#S2tc = np.array([45.8, 46.3, 37.6])

# set up 'true' parameters (these are based on real 0% glycerol measurement)
omega = np.array([[-196.56, 306.62, -708.39]])
amplitude = np.array([2.2, 2.06, 2.40])
amplitude = amplitude / np.mean(amplitude)
R2 = np.array([10.9, 10.4, 8.5])
S2tc = np.array([6.2, 5.8, 4.5])

theta = np.concatenate((amplitude, R2, S2tc)).ravel()

theta_0 = np.array([[1, 1, 1],
                   [30, 30, 30],
                   [15,15,15]]).ravel()


# print all parameters
def print_all():
    print('PARAMETERS\tn_spins\t{:d}'.format(n_spins))
    print('PARAMETERS\tnoise\t{:.2f}'.format(noise))
    print('PARAMETERS\tomega\t{:s}'.format(utils.strrep(omega)))
    print('PARAMETERS\tamp\t{:s}'.format(utils.strrep(amplitude)))
    print('PARAMETERS\tR2\t{:s}'.format(utils.strrep(R2)))
    print('PARAMETERS\tS2tc\t{:s}'.format(utils.strrep(S2tc)))
    for spin in range(n_spins):
        print('THETA_0\t{:d}\tamp\t{:.3f}'.format(spin+1,theta_0[spin]))
    for spin in range(n_spins):
        print('THETA_0\t{:d}\tR2\t{:.3f}'.format(spin+1,theta_0[spin+3]))
    for spin in range(n_spins):
        print('THETA_0\t{:d}\tS2tc\t{:.3f}'.format(spin+1,theta_0[spin+6]))


# routine to 'collect' experimental observation (with noise)
def acquire_point(tau, phi=0.):
    yobs = rs.y(tau, phi, theta, omega) + noise * np.random.standard_normal((1, n_spins))
    if VERBOSE:
        print('ACQUIRE\tTAU\t{:.4f}\tPHASE\t{:.1f}\tYOBS\t{:s}'.format(tau,phi*180/np.pi,utils.strrep(yobs)))
    return yobs


def print_estimate(theta_hat, sigma, label='ESTIMATE'):
    for spin in range(n_spins):
        print('{:s}\t{:d}\tamp\t{:.3f}\t{:.3f}\t{:.3f}'.format(label, spin+1,theta[spin],theta_hat[spin],sigma[spin]))
    for spin in range(n_spins):
        print('{:s}\t{:d}\tR2\t{:.3f}\t{:.3f}\t{:.3f}'.format(label, spin+1,theta[spin+3],theta_hat[spin+3],sigma[spin+3]))
    for spin in range(n_spins):
        print('{:s}\t{:d}\tS2tc\t{:.3f}\t{:.3f}\t{:.3f}'.format(label, spin+1,theta[spin+6],theta_hat[spin+6],sigma[spin+6]))


def print_ensemble_analysis(fit_results):
    mean = np.mean(fit_results,axis=0)
    sd = np.std(fit_results,axis=0)
    print_estimate(mean, sd, 'ENSEMBLE_MEAN_SD')

