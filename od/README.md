# Protocol: optimal sampling of relaxation

## Setup

### Python installation

Required libraries:
* numpy
* scipy
* matplotlib
* numba
* json

### Setting up python scripts

Ensure #! in `analyse_results.py` and `final_results.py` points to the correct python interpreter.


### Setting up bruker bits

Copy files from bruker-bits into python script directory:

```
cp run_adaptive_methyl_ccr.py /opt/topspin/exp/stan/nmr/py/user/
cp -r simplejson /opt/topspin/exp/stan/nmr/py/user/
```

Edit `run_adaptive_methyl_ccr.py` to set appropriate experiment names and directories (at top of script file)

### Setup pulse sequence

Use sequence `adaptive_hsqcphpr.cw`:

acquisition parameters:
* o1p = 0.7 ppm
* cnst21 = (bf hz) for off-resonance pre-sat
* d0 = 13C evolution time (compensating for 13C 90 pulses) > 20 us
* phcor7 = detection phase
* o2p = 11 ppm
* aq = 100 ms, td = 4k
* d1 = 1 s

processing parameters:
* WDW = EM
* LB = 5 Hz
* COROFFS = cnst21 - o1 (i.e. put solvent suppression on-resonance with water)
* BCmod = qfil
* SI = 16k
* set phase correction based on run with d0=20us and phcor7=0

### Chemical shifts

Measure using HMQC:

Peak:		
1H shift:	-0.2573	0.1958	0.4057	0.65367	0.7437
13C shift:	10.4281	11.3917	10.9918	10.784	11.8283

### Optional debugging flags

`od/optimal_design.py`: PLOT = True


## Running

`run_adaptive_methyl_ccr`

* monitor output directory for result files
* don't let anyone use the computer during acquisition!

To interrupt acquisition:
* create a file named 'STOP' in the main experiment directory:
  `touch STOP`
 

## Analysing results

Use script `final_results.py`:
* Pass the experiment directory (containing `adapt_ccr_*` files) as a command line argument
* Set the number of seed points at the top of the script file (for plotting)




