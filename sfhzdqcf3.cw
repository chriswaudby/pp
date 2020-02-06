;SOFAST-H(Z/D)QC
;  Waudby, Ouvry, Davis & Christodoulou (J Biomol NMR, in press)
;  for simultaneous collection of SOFAST-HZQC and SOFAST-HDQC spectra
;
;run as pseudo-3D (td1 = 2)
;process using nmrPipe scripts supplied

prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

"d11=30m"
"d12=20u"
"d13=4u"
"d21=1s/(cnst4*2)"

"in0=inf2"
"td1=2"
"l0=1"
"acqt0=de"

# ifdef ONE_D
"d0=0.1u"
#else
"d0=in0/2-p21*4/3.1415"
# endif /*ONE_D*/

"DELTA1=d21-p16-d16-p39*cnst39"
"DELTA2=p39*cnst39-de-4u"
"DELTA3=DELTA1-p40*0.5"

# ifdef OFFRES_PRESAT
  "TAU=d1-10m-60u-d12*2-d13"
# else
  "TAU=d1-10m"
# endif /*OFFRES_PRESAT*/

"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"

aqseq 312

1 ze 
  d11 pl26:f3
2 10m do:f3

# ifdef OFFRES_PRESAT
  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  TAU cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  30u fq=0:f1
# else
  TAU
# endif /*OFFRES_PRESAT*/

  d12 pl3:f3
  50u UNBLKGRAD

  ; purge Nz
  (p21 ph1):f3
  p16:gp2
  d16

  ; begin main sequence
  (p39:sp23 ph11):f1
  p16:gp1
  d16

  if "l0 %2 == 1"
     {
  ; 1H 180 before d0
  (lalign (DELTA3 p40:sp24 ph1) (DELTA1 p21 ph12 d0 p21 ph1 DELTA1):f3 )
     }
  else
     {
  ; 1H 180 after d0
  (ralign (p40:sp24 ph1 DELTA3 ) (DELTA1 p21 ph12 d0 p21 ph1 DELTA1):f3 )
     }

  DELTA2
  p16:gp1
  d16 pl26:f3
  4u BLKGRAD

  go=2 ph31 cpd3:f3 
  10m do:f3 mc #0 to 2 
     F1QF(ip12)
     F2EA(rp12 & iu0, id0)

exit 
 
ph1=0 
ph2=0 
ph11=0 
ph12=0 2 
ph31=0 2

;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl26: f3 channel - power level for CPD/BB decoupling (low power)
;sp23: f1 channel - shaped pulse 120 degree (Pc9_4_120.1000)
;sp24: f1 channel - shaped pulse 180 degree (Reburp.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (2325 us at 800 MHz)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000              (1700 us at 800  MHz)
;d0 : incremented delay (2D) = in0/2-p21*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)NH
;o1 : recommended to place on H(N) = cnst19
;cnst4: = J(NH)
;cnst19: H(N) chemical shift (offset, in ppm)  [8.2 ppm]
;cnst21: frequency (in Hz) for off-resonance presaturation
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;NS: 2 * n
;DS: 16
;aq: <= 100 msec (with low power cpd, max. 50% duty cycle)
;td1: 2
;td2: number of experiments
;FnMODE[1]: QF
;FnMODE[2]: Echo-AntiEcho
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence
;          use pulse of >= 350 usec

;Options:
;  -DOFFRES_PRESATURATION : off-resonance presaturation during d1 (cnst21, pl9)
;  -DONE_D : for 1D measurement with minimal evolution time (td2=1, td1=2)

;use gradient ratio:	gp 1 : gp 2
;			  11 :    7
;
;for z-only gradients:
;gpz1: 11%
;gpz2:  7%
;
;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
