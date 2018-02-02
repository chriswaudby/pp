#!/usr/bin/env python

import numpy as np

def get_input(prompt, default):
    return input(prompt) or str(default)

N = int(get_input('Number of NUWS dimensions [1]: ', 1))
cos_power = int(get_input('Power of window function, n (cos^n) [2]: ',2))
Nmax = int(get_input('Maximum number of repeats [16]: ', 16))
print('Please enter time domain sizes as REAL points')
if N==1:
    td1 = int(get_input('td (nominal) [64]: ', 64)) // 2
    td2 = 1
else:
    td1 = int(get_input('td (outer loop, nominal) [64]: ', 64)) // 2
    td2 = int(get_input('td (inner loop, nominal) [16]: ', 16)) // 2

# calculate effective td sizes
td1eff=td1
td2eff=td2
for i in range(td1):
    if int(round(Nmax * np.cos(i/td1*np.pi/2)**cos_power))==0:
        td1eff = i
        break        
if N>1:
    for i in range(td2):
        if int(round(Nmax * np.cos(i/td2*np.pi/2)**cos_power))==0:
            td2eff = i
            break        

print()
total_cycles = 0

# now calculate vc lists
if N==1: # 2D
    for i in range(td1eff):
        w1 = np.cos(i/td1*np.pi/2)**cos_power
        c = int(round(Nmax*w1))
        print(c)
        print(c) # second copy for complex points
        total_cycles += 2*c
else: # 3D
    for i in range(td1eff): # outer loop
        w1 = np.cos(i/td1*np.pi/2)**cos_power
        for i2 in range(2): # complex pts
            for j in range(td2eff): # inner loop
                w2 = np.cos(j/td2*np.pi/2)**cos_power
                c = int(round(Nmax * w1 * w2))
                print(c)
                print(c) # second copy for complex points
                total_cycles += 2*c

print()
print('Effective td1 = ' + str(2*td1eff))
if N>1:
    print('Effective td2 = ' + str(2*td2eff))
print('Total cycles = ' + str(total_cycles))

        

