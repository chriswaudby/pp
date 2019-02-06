'''
Adaptive sampling for methyl CCR measurements
running under Jython (urgh)

Chris Waudby
Jul 2017
'''

import math
import sys
import time
import simplejson as json
from os.path import isfile
from subprocess import *

# debugging flag - plot dispersion function d(x) for each iteration?
PLOTDISPERSION = False

# SETUP BITS
analysis_script_path = '/home/waudbyc/pp_950/od'  # path to scripts for numerical analysis in proper python
data_directory = '/home/waudbyc/nmr'
experiment_name = 'charles_ccr_750ile_171218'
BF2 = 238.9908430 



# split up csv list and convert to floats
def split_csv(inputstring):
    tmp = [x.strip() for x in inputstring.split(',')]
    return [float(x) for x in tmp]

# get experiment
proton_template_expt = 2
diffusion_template_expt = 6
hsqc_template_expt = 5
first_expt = 10
dlg_title = 'Adaptive methyl CCR measurement' # default title for dialogs
result = INPUT_DIALOG(title=dlg_title,
    header='Please input the experiment name and IDs of reference experiments.',
    items=['Experiment folder = ',
     '1D template = ',
     'Adaptive HSQC template = ',
     'Diffusion template (zero => do not run) = ',
     'First output experiment = '
    ],
    values=[experiment_name,
     str(proton_template_expt),
     str(hsqc_template_expt),
     str(diffusion_template_expt),
     str(first_expt)
    ])
if result is None:
    MSG('Setup and acquisition aborted!', title=dlg_title)
    EXIT()
else: # unpack results
    experiment_name = result[0]
    proton_template_expt = int(result[1])
    hsqc_template_expt = int(result[2])
    diffusion_template_expt = int(result[3])  # NB change 3 to 4 below when uncommenting!
    first_expt = int(result[4])
if diffusion_template_expt <= 0:
    run_diffusion = False
else:
    run_diffusion = True

# set working directory for temporary files to current experiment
working_directory = "%s/%s" % (data_directory, experiment_name)

# Utility function to change experiment number
def re(expt):    
    dataset = "%s/%i/pdata/1" % (working_directory, expt)
    print("re: %s" % dataset)
    RE_PATH(dataset, show='y')
    time.sleep(0.5) # delay for disk access (prevent crashes?)

# Utility functions for file transfer to/from json (for python interchange)
def loadvar(name):
    # Load parameter from json file in working directory.
    f = open(working_directory + '/adapt_ccr_' + name, 'r')
    x = json.load(f)
    f.close()
    return x

def savevar(x, name):
    # Save parameter x to json file in working directory.
    f = open(working_directory + '/adapt_ccr_' + name, 'w')
    json.dump(x, f)
    f.close()

# get 13C offset, peak positions, initial parameter estimates and seed evolution times
result = INPUT_DIALOG(title=dlg_title,
    header='Please input peak parameters.',
    items=['13C offset (o2p) = ',
        'Peak positions (1H, ppm, comma separated) = ',
        'Peak positions (13C, ppm, comma separated) = ',
        'Relative amplitudes (comma separated) = ',
        '13C R2 (s-1, comma separated) = ',
        'S2tc (ns, comma separated) = ',
        'Evolution times for priming experiments (ms, comma separated, 0.02ms min) = ',
        '13C excitation phase for priming experiments (0 to 359 degrees) = ',
        'Full width of integration regions (ppm) = '],
    values=['11',
        '0.8497, 0.7597, 0.4956, 0.2964, -0.1709',
        '11.7785, 10.7830, 11.0326, 11.4871, 10.5668',
        '1, 1, 1, 1, 1',
        '15, 15, 15, 15, 15',
        '10, 10, 10, 10, 10',
        '0.00002, 0.001, 0.002, 0.005, 0.01, 0.1',
        '0, 0, 0, 0, 0, 0',
#        '0.00002, 0.005, 0.010, 0.040, 0.100',
#        '0, 0, 0, 0, 0',
        '0.03']
    )
if result is None:
    MSG('Setup and acquisition aborted!', title=dlg_title)
    EXIT()
# unpack results
offset_13C = float(result[0])
peak_positions_1H = split_csv(result[1])
peak_positions_13C = split_csv(result[2])
peak_amplitudes = split_csv(result[3])
peak_lambda = split_csv(result[4])
peak_S2tc = split_csv(result[5])
priming_times = split_csv(result[6])
priming_phases = split_csv(result[7])
integration_width = float(result[8])

print("1H peak positions")
print(peak_positions_1H)
print("13C peak positions")
print(peak_positions_13C)

savevar(peak_positions_1H, "peak_positions_1H")
savevar(peak_positions_13C, "peak_positions_13C")
savevar(integration_width, "integration_width")
savevar(priming_times, "priming_times")
savevar(priming_phases, "priming_phases")

# form initial parameter matrix
theta = [peak_amplitudes,
         peak_lambda,
         peak_S2tc]

# save initial parameter estimates to disk
savevar(theta, 'theta0')

# calculate offsets
omega = []
for pp in peak_positions_13C:
    omega.append(2*math.pi*BF2*(pp-offset_13C))
    
# save offsets to disk for later analysis
savevar(omega, 'omega')


# get parameters for interleaving
if run_diffusion:
    result = INPUT_DIALOG(title=dlg_title,
        header='Please input the number of relaxation measurements between interleaving 1D and diffusion experiments.',
        items=['CCR block size = '],
        values=[str(10)])
    if result is None:
        MSG('Setup and acquisition aborted!', title=dlg_title)
        EXIT()
    else: # unpack results
        interleave_block_size = int(result[0])
else:
    interleave_block_size = 0

# get total number of experiments
result = INPUT_DIALOG(title=dlg_title,
    header='Please input the total number of CCR observations to acquire.',
    items=['Number of experiments (incl. seed) = '],
    values=[str(10)])
if result is None:
    MSG('Setup and acquisition aborted!', title=dlg_title)
    EXIT()
N_experiments = int(result[0]) - len(priming_times)

# confirm before starting acquisition
if CONFIRM(message="Begin acqusition?", title=dlg_title) == 0:
    EXIT()

# prepare empty lists for storing experiments and associated evolution times
expts = []
taus = []
phases = []
integrals = []

def prepare_ccr_expt(expt, tau, phase=0):
    print("prepare_ccr_expt: expt=%i, tau=%g, phase=%g" % (expt, tau, phase))
    re(hsqc_template_expt)
    XCMD("wrpa %i" % expt)
    time.sleep(0.5) # delay for disk access (prevent crashes?)
    print("wrpa %i done!" % expt)
    re(expt)
    print("PUTPAR D 0 <%s>" % str(tau))
    PUTPAR('D 0', str(tau))
    print("PUTPAR PHCOR 7 <%s>" % str(phase))
    PUTPAR('PHCOR 7', str(phase))
    print("PUTPAR O2P <%s>" % str(offset_13C))
    PUTPAR('O2P', str(offset_13C))
    expts.append(expt)
    taus.append(tau)
    phases.append(phase)
    print("experiment prepared!")

def prepare_diffusion_expt(expt):
    print("prepare_diffusion_expt: expt=%i" % (expt))
    re(diffusion_template_expt)
    XCMD("wrpa %i" % expt)
    time.sleep(0.5) # delay for disk access (prevent crashes?)
    print("wrpa %i done!" % expt)
    re(expt)

def prepare_proton_expt(expt):
    print("prepare_proton_expt: expt=%i" % (expt))
    re(proton_template_expt)
    XCMD("wrpa %i" % expt)
    time.sleep(0.5) # delay for disk access (prevent crashes?)
    print("wrpa %i done!" % expt)
    re(expt)

def analyse_expt(expt):
    re(expt)
    EFP()
    # TODO - baseline correction
    y = []
    for x in peak_positions_1H:
        from_ppm = x - 0.5*integration_width
        to_ppm = x + 0.5*integration_width
        spec = GETPROCDATA(from_ppm, to_ppm)
        y.append(sum(spec))
    integrals.append(y)
    print('Peak integrals:')
    print(y)

# stop acquisition if a file 'STOP' is detected in working directory
def stop():
    if isfile(working_directory + '/STOP'):
        return True
    else:
        return False

current_expt = first_expt - 1
total_relaxation_expts = 0 # counter for interleaving diffusion

for tau, phase in zip(priming_times, priming_phases):
    if stop():
        print('STOP during seed experiments, current_expt = %i' % current_expt)
        break
    current_expt += 1
    prepare_ccr_expt(current_expt, tau, phase)
    print('Running seed experiment %i: tau = %g, phase = %g' % (current_expt, tau, phase))
    #XCMD('zg', WAIT_TILL_DONE) # BUG - this doesn't actually wait until execution is finished!
    XCMD('au_zgonly', WAIT_TILL_DONE) # solution - run indirectly via AU program
    print('Finished seed experiment %i!' % (current_expt))
    analyse_expt(current_expt)
    total_relaxation_expts += 1

for iteration in range(N_experiments):
    if stop():
        print('STOP during main experiments, current_expt = %i' % current_expt)
        break

    # write current status to disk for analysis in proper python
    savevar(integrals, 'integrals')
    savevar(taus, 'taus')
    savevar(phases, 'phases')
    savevar(theta, 'theta')
    
    # calculate theta update and next sampling point
    output = Popen(["%s/analyse_results.py" % analysis_script_path,
                    working_directory], \
                   stdout=PIPE).communicate()[0]
    print("Fit results:")
    print(output)
    
    # load analysis results
    tau = loadvar('newtau')
    phase = loadvar('newphase')
    print("Optimised control point: tau = %s, phase = %s" % (tau,phase))
    theta = loadvar('newtheta')
    print("Updated parameter vector:")
    print(theta)
    
    # show dispersion curve
    #if PLOTDISPERSION:
    #    dispersioncurve = loadvar('dispersioncurve')
    #    dispersioncurve_t, dispersioncurve_m = zip(*dispersioncurve)
    #    props = GET_DISPLAY_PROPS(-0.005, 0.155, "s","Evolution time","Dispersion","Adaptive methyl CCR sampling: dispersion curve","line")
    #    DISPLAY_DATALIST_XY([dispersioncurve_m],[dispersioncurve_t],[props],"Methyl CCR measurements")
    
    # now acquire the next and optimal control point
    current_expt += 1
    prepare_ccr_expt(current_expt, tau, phase)
    print('Running experiment %i: tau = %g, phase = %g' % (current_expt, tau, phase))
    XCMD('au_zgonly', WAIT_TILL_DONE)
    print('Finished experiment %i!' % (current_expt))
    analyse_expt(current_expt)

    # increment expt counter for interleaving
    total_relaxation_expts += 1

    if run_diffusion and (total_relaxation_expts % interleave_block_size == 0):
        current_expt += 1
        prepare_proton_expt(current_expt)
        print('Running proton 1D, experiment %i' % (current_expt))
        XCMD('au_zgonly', WAIT_TILL_DONE)
        print('Finished 1D experiment %i!' % (current_expt))
        current_expt += 1
        prepare_diffusion_expt(current_expt)
        print('Running diffusion, experiment %i' % (current_expt))
        XCMD('au_zgonly', WAIT_TILL_DONE)
        print('Finished diffusion experiment %i!' % (current_expt))

# write current status to disk for analysis in proper python
savevar(integrals, 'integrals')
savevar(taus, 'taus')
savevar(phases, 'phases')
savevar(theta, 'theta')

# calculate theta update and next sampling point
output = Popen(["%s/analyse_results.py" % analysis_script_path,
                working_directory], \
               stdout=PIPE).communicate()[0]
print("Fit results (final):")
print(output)

print("FINISHED!")

