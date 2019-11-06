; SOFAST-HMQC for measurement of methyl MQ relaxation
; using purge element to suppress outer lines
;
; Removal of 13C equilibrium magnetisation (for methyl TROSY)
; Addition of clean-up gradient-pair
; Delays adjusted for zero first-order phase correction


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define delay XI1
define delay XI2
define delay vdMin

"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d11=30m"
"d12=20u"
"d13=4u"

"TAU=d1-d11-4u-d12-50u-p3-d13-p16-d16"
"DELTA1=d2-p39*cnst39-p3"
"XI1=d2/4-p17-d17-0.5*p40"
"XI2=d2/4-p17-d17-p3-0.5*p40"
"DELTA2=d2-d12-4u-p40*0.5"
"vdMin=larger(4*p19+4*d17,2*d0+2*p40)+4*p3"

"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"

"in0=inf2"
"d0=in0/2"

"acqt0=0"
baseopt_echo
aqseq 312


1 ze 
  vdMin
  d11 pl12:f2
2 d11 do:f2
  4u BLKGRAD


  ; recycle delay - calculate new delays for relaxation in this time
  TAU
  "d20=vd/4-p3-p19-d17"
  "d21=vd/4-p3"
  "d22=vd/4-p3+d0*0.5-0.5*p40"
  "d23=vd/4-p3-p40*0.5"

  d12 pl1:f1 pl2:f2
  50u UNBLKGRAD

  ; crush eq'm 13C magnetisation
  (p3 ph1):f2    
  d13
  p16:gp0
  d16

  (p39:sp23 ph1):f1
  ;"DELTA1=d2-p39*cnst39-p3"
  DELTA1

  (p3 ph11):f2

  ;begin purge element
  ;"XI1=d2/4-p17-d17-0.5*p40"
  XI1
  p17:gp1
  d17
  (center (p40:sp24 ph1):f1 (p4 ph1):f2 )
  p17:gp1
  d17
  ;"XI2=d2/4-p17-d17-p3-0.5*p40"
  XI2
  (p3 ph12):f2

  ;begin relaxation/t1 evolution
  ;"d20=vd/4-p3-p19-d17"
  d20
  p19:gp2
  d17
  (p4 ph1):f2
  ;"d21=vd/4-p3"
  ;"d22=vd/4-p3+d0*0.5-0.5*p40"
  (lalign (d22 p40:sp24 ph1):f1 (d21 d0 d20):f3 )
  p19:gp2
  d17
  (p4 ph1):f2
  ;"d23=vd/4-p3-p40*0.5"
  d23
  (p40:sp24 ph13):f1 (p3 ph14):f2

  ;back-transfer
  ;"DELTA2=d2-d12-4u-p40*0.5"
  DELTA2
  d12 pl12:f2
  4u BLKGRAD

  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2
       F1QF(ivd)
       F2PH(ip14, id0)
  4u BLKGRAD
exit 
  

ph1= 0 
ph11=0 2
ph12=1 1 3 3
ph13=0 0 0 0 2 2 2 2
ph14=0
ph29=0
ph31=0 2

;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;sp23: f1 channel - shaped pulse 120 degree
;                   (Pc9_4_120.1000 or Q5.1000)
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse
;p17: gradient pulse [300 usec]
;p19: gradient pulse [200 usec]
;p22 : f3 channel - 180 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (3.0ms at 600.13 MHz)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000               (1.0ms at 600.13 MHz)
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery [200 usec]
;d17: short delay for gradient recovery [100 usec]
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;cnst19: H(N) chemical shift (offset, in ppm)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;           Q5.1000: -0.07
;nd1: 1
;NS: 8 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz0: -40%
;gpz1: 31%
;gpz2: 11%

;use gradient files:
;gpnam0: SMSQ10.100
;gpnam1: SMSQ10.32
;gpnam2: SINE.20

                                          ;preprocessor-flags-start
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
;	    option -DOFFRES_PRESAT (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


