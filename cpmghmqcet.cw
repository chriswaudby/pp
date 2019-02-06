;CPMG-HMQC with gradient selection
;hmqcet.cw
;1H,13C HMQC
;phase sensitive using Echo/Anti-Echo (Hurd & John 1991)
; Modern Instrumental Analysis, edited by Satinder Ahuja, Neil Jespersen, p286
;with decoupling during acquisition
;
;modified Chris Waudby 25/11/18


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

"d11=30m"
"d12=20u"
"d21=1s/(cnst2*2)"
"p2=p1*2"

"in0=inf1/2"

"d0=3u"
"DELTA1=d21/32-larger(p1,p3)"
"DELTA=DELTA1*2"

"DELTA2=p17+d17-p1*0.6366"
"TAU=p3*0.6366+p17+d17-p1-d0"
"acqt0=0"
baseopt_echo

1 ze
  d11 pl12:f2
2 d1 do:f2
  d12 pl1:f1 pl2:f2
  50u UNBLKGRAD

  ; purge equilibrium 13C
  (p3 ph1):f2
  4u
  p16:gp3
  d16

  ; begin main sequence
  (p1 ph1):f1

#ifdef CPMG1
; phases: x y x y  y x y x  -x -y -x -y  -y -x -y -x
;         1 2 1 2  2 1 2 1  3  4  3  4   4  3  4  3
DELTA1
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA1
#else
  d21
#endif /* CPMG1 */

  ; begin t1 period
  (p3 ph13):f2

  ; first spin-echo for gradient selection
  p17:gp1
  d17
  (p4 ph1):f2
  TAU

  ; t1 and central refocusing pulse
  d0
  (p2 ph12):f1
  d0

  ; second spin-echo for gradient selection
  TAU
  (p4 ph1):f2
  p17:gp1
  d17
  (p3 ph14):f2
  ; end t1 period

#ifdef CPMG2
; phases: x y x y  y x y x  -x -y -x -y  -y -x -y -x
;         1 2 1 2  2 1 2 1  3  4  3  4   4  3  4  3
DELTA1
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA1
#else
  d21
#endif /* CPMG2 */

  ; extra spin-echo for gradient selection
  p17:gp2*EA
  d17 pl12:f2
  (p2 ph2):f1
  DELTA2 BLKGRAD

  ; acqusition
  go=2 ph31 cpd2:f2
  d1 do:f2 mc #0 to 2
     F1EA(igrad EA, id0 & ip13*2 & ip31*2)
exit


ph1=0
ph2=1
ph3=2
ph4=3
ph12=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph13=0 2
ph14=0 0 2 2
ph31=0 2 2 0 2 0 0 2


;pl2 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling (low power)
;p16: homospoil/gradient pulse                       [1 msec]
;p1:  f1 channel -  90 degree high power pulse
;p3:  f2 channel -  90 degree high power pulse
;d0 : incremented delay (2D) = in0/2-p3*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)NH
;cnst2: = J(CH) [125 Hz for methyls]
;inf1: 1/SW(C) = 2 * DW(C)
;in0: 1/ SW(C) = 2 * DW(C)
;nd0: 1
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: Echo-AntoEcho
;cpd2: decoupling according to sequence defined by cpdprg2: garp4
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:	gp 1 : gp 2
;			  40 :    20.1


;for z-only gradients:
;gpz1: 40%
;gpz2: 20.1%
;gpz3: 31%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100

