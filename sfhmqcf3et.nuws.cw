;sfhmqcf3et.nuws.cw
;1H,15N SOFAST HMQC
;phase sensitive using Echo/Anti-Echo (Hurd & John 1991)
; Modern Instrumental Analysis, edited by Satinder Ahuja, Neil Jespersen, p286
;with apodisation-weighted sampling
;with decoupling during acquisition
;
;set sampling schedule via vclist:
;  total number of scans = ns * c
;  add one to first point of vclist to allow for dummy scans
;
;modified Chris Waudby 02/05/18
;
;P.Schanda and B. Brutscher, J. Am. Chem. Soc. 127, 8014 (2005)
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

#ifdef NUWS
define loopcounter myCounter
"myCounter=0"
#endif /* NUWS */

"cnst4=92"

"d11=30m"
"d12=20u"
"d21=1s/(cnst4*2)"

"p22=2*p21"

"in0=inf1/2"

"d0=3u"

"DELTA=p21*0.6366+p17+d17-p40*0.5-d0"
"DELTA1=d21-p39*cnst39"
"DELTA2=d21-p17-d16-4u"
"acqt0=0"


"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"


1 ze
  d11 pl26:f3
2 d1 do:f3
3 d12 pl3:f3
  50u UNBLKGRAD

  ; purge equilibrium 15N
  (p21 ph1):f3
  4u
  p16:gp3
  d16

  (p39:sp23 ph1):f1
  DELTA1

  (p21 ph3):f3
  p17:gp1
  d17
  (p22 ph2):f3
  DELTA

  d0
  (p40:sp24 ph1):f1
  d0

  DELTA
  (p22 ph1):f3
  p17:gp1
  d17
  (p21 ph4):f3

  DELTA2
  p17:gp2*EA
  d16 pl26:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3

#ifdef NUWS
  if "myCounter==1" goto 10
  zd
  "myCounter=1"
10 4u
  ; repeat acquisition block according to schedule in vclist
  lo to 2 times c
  3m ivc
#endif /* NUWS */

  d1 do:f3 mc #0 to 2
     F1EA(igrad EA, id0 & ip3*2 & ip31*2)
exit


ph1=0
ph2=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph3=0 2
ph4=0 0 2 2
ph31=0 2 2 0 2 0 0 2


;pl3 : f3 channel - power level for pulse (default)
;pl26: f3 channel - power level for CPD/BB decoupling (low power)
;sp23: f1 channel - shaped pulse 120 degree
;                   (Pc9_4_120.1000 or Q5.1000)
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p21: f3 channel -  90 degree high power pulse
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
;cnst4: = J(NH) [92 Hz]
;cnst19: H(N) chemical shift (offset, in ppm)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/ SW(N) = 2 * DW(N)
;nd0: 1
;NS: 4 * n
;DS: 16
;aq: <= 50 msec
;td1: number of experiments
;FnMODE: Echo-AntiEcho
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence
;          use pulse of >= 350 usec

;for z-only gradients:
;gpz1: 40%
;gpz2: -8.1%
;gpz3: 31%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100


                                          ;preprocessor-flags-start
                                          ;preprocessor-flags-end


;Processing

;PHC0(F1): 0
;PHC1(F1): 0
;FCOR(F1): 0.5

