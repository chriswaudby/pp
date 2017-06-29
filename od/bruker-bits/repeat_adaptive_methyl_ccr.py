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

N_repeats = 10

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
     'First output experiment = '
    ],
    values=[experiment_name,
     str(proton_template_expt),
     str(hmqc_template_expt),
     str(first_expt)
    ])
if result is None:
    MSG('Setup and acquisition aborted!', title=dlg_title)
    EXIT()
else: # unpack results
    experiment_name = result[0]
    proton_template_expt = int(result[1])
    hmqc_template_expt = int(result[2])
    first_expt = int(result[3])

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
        '15, 15, 15, 15, 15',
        '10, 10, 10, 10, 15',
#        '0.00002, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04',
#        '0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0',
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

################################################################################
# BEGIN REPEATS
################################################################################
for repeat in range(N_repeats):

    # set working directory for temporary files to current experiment and repeat
    working_directory = "/opt/topspin3.5pl2/data/jc/nmr/%s_%i/" % (experiment_name, repeat)

    ## form initial parameter matrix
    theta = [peak_amplitudes,
             peak_lambda,
             peak_S2tc]


    # Utility function to change experiment number
    def re(expt):    
        dataset = "/opt/topspin3.5pl2/data/jc/nmr/%s_%i/%i/pdata/1" % (experiment_name, repeat, expt)
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

    # save initial parameter estimates to disk
    savevar(theta, 'theta0')

    # calculate offsets
    omega = []
    for pp in peak_positions_13C:
        omega.append(2*math.pi*BF2*(pp-offset_13C))
        
    # save offsets to disk for later analysis
    savevar(omega, 'omega')


    # prepare empty lists for storing experiments and associated evolution times
    expts = []
    taus = []
    phases = []
    integrals = []

    def prepare_ccr_expt(expt, tau, phase=0):
        re(hmqc_template_expt)
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

