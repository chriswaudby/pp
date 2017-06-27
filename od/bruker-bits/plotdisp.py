'''
Adaptive sampling for methyl CCR measurements

Chris Waudby
Jan 2016
'''

import math
import sys
import simplejson as json
from subprocess import *
#import de.bruker.nmr.prsc.toplib as top

BF2 = 176.083032

dlg_title = 'Adaptive methyl CCR measurement' # default title for dialogs

# split up csv list and convert to floats
def split_csv(inputstring):
    tmp = [x.strip() for x in inputstring.split(',')]
    return [float(x) for x in tmp]

# get experiment
experiment_name = 'chris_adaptive_methyl_ccr_090116'
result = INPUT_DIALOG(title=dlg_title,
    header='Please input the experiment name and IDs of reference experiments.',
    items=['Experiment folder = '],
    values=[experiment_name])
if result is None:
    MSG('Setup and acquisition aborted!', title=dlg_title)
    EXIT()

experiment_name = result[0]

working_directory = "/opt/topspin3.2/data/jc/nmr/%s/" % experiment_name
def loadvar(name):
    # Load parameter from json file in working directory.
    f = open(working_directory + 'adapt_ccr_' + name, 'r')
    x = json.load(f)
    f.close()
    return x


# show dispersion curve
dispersioncurve = loadvar('dispersioncurve')
dispersioncurve_t, dispersioncurve_m = zip(*dispersioncurve)
#print dispersioncurve_t
#print dispersioncurve_m
props = GET_DISPLAY_PROPS(-0.005, 0.155, "s","Evolution time","Dispersion","Adaptive methyl CCR sampling: dispersion curve","line")
DISPLAY_DATALIST_XY([dispersioncurve_m],[dispersioncurve_t],[props],"Methyl CCR measurements")

