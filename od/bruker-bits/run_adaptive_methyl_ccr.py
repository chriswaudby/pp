'''
Adaptive sampling for methyl CCR measurements
running under Jython (urgh)

Chris Waudby
Jul 2017
'''

import math
import sys
import simplejson as json
from subprocess import *

# debugging flag - plot dispersion function d(x) for each iteration?
PLOTDISPERSION = False

# path to scripts for numerical analysis in proper python
analysis_script_path = '/home/nmrsu/pp_cw/od/'

BF2 = 201.2230230

dlg_title = 'Adaptive methyl CCR measurement' # default title for dialogs

# split up csv list and convert to floats
def split_csv(inputstring):
    tmp = [x.strip() for x in inputstring.split(',')]
    return [float(x) for x in tmp]

# get experiment
experiment_name = 'chris_ccr_750ile_adapt_270717'
proton_template_expt = 2
diffusion_template_expt = 999
hmqc_template_expt = 7
first_expt = 20
result = INPUT_DIALOG(title=dlg_title,
    header='Please input the experiment name and IDs of reference experiments.',
    items=['Experiment folder = ',
     '1D template = ',
     'HMQC template = ',
#     'Diffusion template = ',
     'First output experiment = '
    ],
    values=[experiment_name,
     str(proton_template_expt),
     str(hmqc_template_expt),
#     str(diffusion_template_expt),
     str(first_expt)
    ])
if result is None:
    MSG('Setup and acquisition aborted!', title=dlg_title)
    EXIT()
else: # unpack results
    experiment_name = result[0]
    proton_template_expt = int(result[1])
    hmqc_template_expt = int(result[2])
#    diffusion_template_expt = int(result[3])  # NB change 3 to 4 below when uncommenting!
    first_expt = int(result[3])

# set working directory for temporary files to current experiment
working_directory = "/opt/topspin3.5pl2/data/jc/nmr/%s/" % experiment_name

# Utility function to change experiment number
def re(expt):    
    dataset = "/opt/topspin3.5pl2/data/jc/nmr/%s/%i/pdata/1" % (experiment_name, expt)
    print("re: %s" % dataset)
    RE_PATH(dataset, show='y')

# Utility functions for file transfer to/from json (for python interchange)
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
        'Width of integration regions (ppm) = '],
    values=['11',
        '0.7208, 0.6319, 0.3828, 0.1734, -0.2803',
        '11.8041, 10.76, 10.967, 11.3704, 10.4062',
        '1, 1, 1, 1, 1',
        '15, 15, 15, 15',
        '10, 10, 10, 10',
        '0.00002, 0.005, 0.010, 0.040, 0.100',
        '0, 0, 0, 0, 0',
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
#result = INPUT_DIALOG(title=dlg_title,
#    header='Please input the time interval for interleaving 1D and diffusion experiments.',
#    items=['CCR block duration (minutes) = '],
#    values=[str(90)])
#if result is None:
#    MSG('Setup and acquisition aborted!', title=dlg_title)
#    EXIT()
#else: # unpack results
#    interleave_time = float(result[0])

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
    re(hmqc_template_expt)
    #dataset = "/opt/topspin3.2/data/jc/nmr/%s/%i/pdata/1" % (experiment_name, expt)
    #print("Writing dataset: %s" % dataset)
    #print(type(dataset))
    #dataset=str(dataset)
    #print(type(dataset))
    #print(dataset)
    #top.Cmd.wr(dataset, "n")
    ##WR_PATH(dataset) # , confirm="n") #BUG - override not a valid option?

    XCMD("wrpa %i" % expt)
    print("wrpa %i done!" % expt)
    re(expt)
    PUTPAR('D 0', str(tau))
    PUTPAR('PHCOR 7', str(phase))
    PUTPAR('O2P', str(offset_13C))
    expts.append(expt)
    taus.append(tau)
    phases.append(phase)

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


current_expt = first_expt - 1

for tau, phase in zip(priming_times, priming_phases):
    current_expt += 1
    prepare_ccr_expt(current_expt, tau, phase)
    print('Running seed experiment %i: tau = %g, phase = %g' % (current_expt, tau, phase))
    #XCMD('zg', WAIT_TILL_DONE) # BUG - this doesn't actually wait until execution is finished!
    XCMD('au_zgonly', WAIT_TILL_DONE) # solution - run indirectly via AU program
    print('Finished seed experiment %i!' % (current_expt))
    analyse_expt(current_expt)

for iteration in range(N_experiments):
    # write current status to disk for analysis in proper python
    savevar(integrals, 'integrals')
    savevar(taus, 'taus')
    savevar(phases, 'phases')
    savevar(theta, 'theta')
    
    # calculate theta update and next sampling point
    output = Popen(["/home/nmrsu/pp_cw/od/analyse_results.py",
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
    print('Running experiment %i: tau = %g' % (current_expt, tau))
    XCMD('au_zgonly', WAIT_TILL_DONE)
    print('Finished experiment %i!' % (current_expt))
    analyse_expt(current_expt)

# write current status to disk for analysis in proper python
savevar(integrals, 'integrals')
savevar(taus, 'taus')
savevar(phases, 'phases')
savevar(theta, 'theta')

# calculate theta update and next sampling point
output = Popen(["/home/nmrsu/pp_cw/od/analyse_results.py",
                working_directory], \
               stdout=PIPE).communicate()[0]
print("Fit results (final):")
print(output)

print("FINISHED!")

