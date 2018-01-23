import numpy as np

N = int(input('Number of points in difframp [8]: ') or '8')
G0 = float(input('First gradient point [5%]: ') or '5')
G1 = float(input('Last gradient point [95%]: ') or '95')
delta = float(input('Little delta (ms) [4 ms]: ') or '4') * 0.001
DELTA = float(input('Big delta (ms) [100 ms]: ') or '100') * 0.001
Gmax = float(input('Gmax (T/m) [0.55 T/m]: ') or '0.55')
T = float(input('Temperature (K) [298 K]: ') or '298')
solvent = input('(H)2O or (D)2O [H]: ') or 'H'
if solvent=='D':
    h2o = False
else:
    h2o = True

gammaH = 2.675e8
sigma = 0.9 # SMSQ
k = 1.38e-23

q = (gammaH * delta * sigma * Gmax)**2 * (DELTA - delta/3.)
G = np.linspace(G0*0.01,G1*0.01,N)
G2 = G**2
q = q * G2

output = 'dosyView.tcl -tau '
for val in q:
    output += str(val) + ' '

print()
print('DOSYVIEW COMMAND:')
print(output)

# viscosities: Cho et al, J Phys Chem B (1999) 103 1991-1994
if h2o:
    A = 802.25336
    a = 3.4741e-3
    b = -1.7413e-5
    c = 2.7719e-8
    gamma = 1.53026
    T0 = 225.334
else:
    A = 885.60402
    a = 2.799e-3
    b = -1.6342e-5
    c = 2.9067e-8
    gamma = 1.55255
    T0 = 231.832

DT = T - T0
eta = A * (DT + a*DT**2 + b*DT**3 + c*DT**4)**(-gamma)

print()
print('CONVERSION TO HYDRODYNAMIC RADIUS:')
print('Viscosity (mPa s): '+str(eta))
x = k*T / (6*np.pi*eta*0.001) * 1e9
print()
print('rH [nm] = {0:6.3e} / D'.format(x))
print()

