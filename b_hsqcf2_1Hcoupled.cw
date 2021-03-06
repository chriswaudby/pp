;b_hsqcf2_1Hcoupled.cw
;1H,13C BEST-HSQC with 1H coupling
;phase sensitive
;evolution times adjusted for zero 1H phase correction (AQMOD digital)
; and 1/2 dwell (90/180) 13C phase correction
;with decoupling during acquisition
;
;modified Chris Waudby 15/10/16
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"d11=30m"
"d12=20u"
"d21=1s/(cnst2*4)"

"p4=2*p3"

"in0=inf1"
"d0=in0/2-p3*4/3.1415"


"DELTA1=d21-p17-d16-0.5*larger(p40,p4)-p39*cnst39"
"DELTA2=d21-p17-d16-0.5*larger(p40,p4)-p39*cnst39"
"DELTA3=d21-p17-d16-0.5*larger(p40,p4)-de-8u"


"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"


1 ze
  d11 pl12:f2
2 d1 do:f2
3 d12 pl2:f2
  50u UNBLKGRAD

  ; purge equilibrium 13C
  (p3 ph1):f2
  4u
  p16:gp5
  d16

  (p39:sp23 ph1):f1
  p17:gp1
  d16
  DELTA1
  ;(center (p40:sp24 ph1):f1 (p3 ph11 p4 ph12 p3 ph11):f2 ) ; composite 13C inversion
  (center (p40:sp24 ph1):f1 (p3*2 ph1):f2 )
  DELTA1
  p17:gp1
  d16
  (p39:sp23 ph2):f1

  ; zz purge
  p16:gp3
  200u

  ; t1 evolution
  (p3 ph3):f2
  d0
  (p3 ph4):f2

  ; zz purge
  p16:gp4
  200u

  (p39:sp23 ph1):f1
  p17:gp2
  d16
  DELTA2
  (center (p40:sp24 ph1):f1 (p3*2 ph1):f2 )
  DELTA3
  p17:gp2
  d16
  4u BLKGRAD
  4u pl12:f2
  go=2 ph31 cpd2:f2
  d1 do:f2 mc #0 to 2
     F1PH(ip3, id0)
exit


ph1=0
ph2=1
ph3=0 2
ph4=0 0 2 2
ph11=0
ph12=1
ph31=0 2 2 0


;pl2 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling (low power)
;sp23: f1 channel - shaped pulse 120 degree
;                   (Pc9_4_120.1000 or Q5.1000)
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p3:  f2 channel -  90 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (3.0ms at 600.13 MHz)
;                  (or Q5.1000 (90o)            (2.0ms at 600.13 MHz) )
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000               (1.0ms at 600.13 MHz)
;d0 : incremented delay (2D) = in0/2-p3*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)NH
;cnst2: = J(CH) [125 Hz for methyls]
;cnst19: H(N) chemical shift (offset, in ppm)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;           Q5.1000: -0.07
;inf1: 1/SW(C) = 2 * DW(C)
;in0: 1/ SW(C) = 2 * DW(C)
;nd0: 1
;NS: 4 * n
;DS: 16
;aq: <= 50 msec
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEC
;cpd2: decoupling according to sequence defined by cpdprg2: garp4
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence



;for z-only gradients:
;gpz1: 11%
;gpz2:  7%
;gpz3:  9%
;gpz4: 17%
;gpz5: 13%

;use gradient files:
;gpnam1: SINE.20
;gpnam2: SINE.20
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100


                                          ;preprocessor-flags-start
                                          ;preprocessor-flags-end


;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1



;$Id: sfhmqcf3gpph,v 1.1.2.8 2009/11/18 11:19:58 ber Exp $
