;methyl-SOFAST-HMQC with L2 filter
;Chris Waudby April 2020

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf1"
"in20=inf1*0.5"
"d0=in0/2"

define delay XI1
define delay XI2
"TAU=d1-d11-20u-d12-50u-p3-d13-p16-d16"
"DELTA1=d2-cnst39*p39-0.5*p40-p16-d16"
"DELTA2=d2-p16-d16-d12-4u-de"
"XI1=d2*0.25-p19-d19-p40"
"XI2=d2*0.25-p19-d19-0.5*p40-p3"
"d20=d2*0.25-p19-d19-p40+0.5*d0"

"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"

"acqt0=de"
baseopt_echo

1 ze 
  vdmin
  d11 pl12:f2
2 d11 do:f2
  ; relaxation period
  ;"TAU=d1-d11-d12-50u-p3-d13-p16-d16"
  TAU 
  d12 pl1:f1 pl2:f2
  50u UNBLKGRAD

  (p3 ph1):f2    ; crush eq'm 13C magnetisation
  d13
  p16:gp1
  d16

  ; start main sequence
  (p39:sp23 ph10):f1
  p16:gp2
  d16
  DELTA1

  (center (p40:sp24 ph1):f1 (p3 ph11):f2 )
  p19:gp3
  d19
  XI1
  (center (p40:sp24 ph1):f1 (p4 ph1):f2 )
  p19:gp3
  d19
  (lalign (d20 p40:sp24 ph1):f1 (XI2 p3 ph12 d0 p3 ph13 DELTA2):f2 )
  p16:gp2
  d16
  d12 pl12:f2
  4u BLKGRAD

  ; acquisition
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2
       F1PH(ip13, id0 & id20)

  4u do:f2
  4u BLKGRAD
exit 
  
ph1= 0
ph2= 1 
ph10=0
ph11=0 2
ph12=1 1 3 3
ph13=0 0 0 0 2 2 2 2
ph29=0
ph31=0 2 0 2 2 0 2 0

;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;sp23: f1 channel - shaped pulse 120 degree
;                   (Pc9_4_120.1000)
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse
;p19: gradient pulse [100 usec]
;p22 : f3 channel - 180 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (1328 usec at 800 MHz)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000               (669 usec at 800 MHz)
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery          [200 usec]
;d19: short delay for gradient recovery [200 usec]
;cnst2: = J(CH)
;cnst19: H(Me) chemical shift (offset, in ppm) [0.5 ppm]
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;l0: number of repeats for entire experiment
;NS: 8 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz2: 11%
;gpz3: 40%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.32
