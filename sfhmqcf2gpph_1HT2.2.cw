;SOFAST-HMQC for methyl 1H T2 measurement
;with L2 filter

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf2"
"d0=in0/2-1.2732*p3"

define delay XI1
define delay XI2
"TAU=d1-d11-20u-d12-50u-p3-d13-p16-d16"
"DELTA1=d2-cnst39*p39"
"DELTA2=d2-d12-4u-de"
"XI1=d2*0.25-p40*0.5-p19-d19"
"XI2=d2*0.25-p3-p40*0.5-p19-d19"

define delay vdmin
"vdmin=4*(p1+p3*2+4u+p17+d17)+2*p40"

"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"

"acqt0=de"
baseopt_echo

aqseq 312

1 ze 
  vdmin
  d11 pl12:f2
2 d11 do:f2

  20u
  "d3=d2*0.25-p40-p19-d19+d0*0.5"
  "d20=vd/4-p3*2"
  "d21=vd/4-p3*2-p17-d17-0.5*p40"
  "d22=vd/4-p3*2-p17-d17-0.5*p40"
  "d23=vd/4-p3*2"

  ; relaxation period
  ;"TAU=d1-d11-20u-d12-50u-p3-d13-p16-d16"
  TAU 
  d12 pl1:f1 pl2:f2
  50u UNBLKGRAD

  (p3 ph1):f2    ; crush eq'm 13C magnetisation
  d13
  p16:gp1
  d16

  ; start main sequence
  (p39:sp23 ph10):f1  ; INEPT
  ;"DELTA1=d2-p39*cnst39"
  DELTA1

  (p3 ph11):f2
  p19:gp3
  d19
  ;"XI1=d2*0.25-p40*0.5-p19-d19"
  XI1
  (center (p40:sp24 ph1):f1 (p3 ph1 p4 ph2 p3 ph1):f2 )
  p19:gp3
  d19
  ;"XI2=d2*0.25-p3-p40*0.5-p19-d19"
  ;"d3=d2*0.25-p40-p19-d19+d0*0.5"
  ;"d20=vd/4-p3*2"
  (lalign (d3 p40:sp24 ph14):f1 (XI2 p3 ph12 d0 p3 ph13 d20):f2 )
  (p3 ph1):f2
  (p4 ph2):f2
  (p3 ph1):f2
  ;"d21=vd/4-p3*2-0.5*p40-p17-d17"
  d21
  p17:gp4
  d17
  (p40:sp24 ph1):f1
  p17:gp4
  d17
  ;"d22=vd/4-p3*2-p17-d17-0.5*p40"
  d22
  (p3 ph1):f2
  (p4 ph2):f2
  (p3 ph1):f2
  ;"d23=vd/4-p3*2"
  d23

  ; back-transfer
  d12 pl12:f2
  4u BLKGRAD
  DELTA2

  ; acquisition
  go=2 ph31 cpd2:f2 

  d11 do:f2 mc #0 to 2
       F1QF(ivd)
       F2PH(ip13, id0)

  ; repeat whole experiment
  lo to 2 times l0
  4u do:f2
  4u BLKGRAD
exit 
  
ph1= 0
ph2= 1 
ph10=0
ph11=0 2
ph12=1 1 3 3
ph13=0 0 0 0 2 2 2 2
ph14=(3) {{0}*8}^1^2
ph29=0
ph31=(6) {{{0 3 0 3}^2}^2^4}


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
;p17: gradient pulse [100 usec]
;p19: gradient pulse [50 usec]
;p22 : f3 channel - 180 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (3.0ms at 600.13 MHz)
;                  (or Q5.1000 (90o)            (2.0ms at 600.13 MHz) )
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000               (1.0ms at 600.13 MHz)
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d3 : 1/(8J)CH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;d17: shorter delay for gradient recovery [100 usec]
;d19: short delay for gradient recovery [100 usec]
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;cnst19: H(Me) chemical shift (offset, in ppm)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;l0: number of repeats for entire experiment
;NS: 8 * n, 24 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz3: -40%
;gpz4: 11%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SINE.10
;gpnam4: SINE.10
                                          ;preprocessor-flags-start
;SINGLEDWELL: for initial sampling delay of one dwell-time with 
;	    option -DSINGLEDWELL (eda: ZGOPTNS)
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
;	    option -DOFFRES_PRESAT (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: hmqcphpr,v 1.4 2007/04/11 13:34:30 ber Exp $
