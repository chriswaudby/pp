;hmqcetgppr.2
;avance-version (12/01/11)
;HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive using Echo/Antiecho gradient selection
;with decoupling during acquisition
;using shaped pulses for inversion and refocussing on f2 - channel
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"d2=1s/(cnst2*2)"
"d11=30m"


"d0=3u"

"in0=inf1/2"


"DELTA=p2+d0*2"
"DELTA1=d2-p16-8u"
"DELTA2=p16+d16+d0+p2/2+p24/2-p14/2-4u"


1 ze 
2 d11 do:f2
 
  4u pl9:f1
  d1 cw:f1 ph1
  4u do:f1
  4u pl1:f1

  p1 ph1
  d2
  DELTA2 pl0:f2
  4u
  (p14:sp3 ph5):f2
  4u
  DELTA2 pl2:f2 UNBLKGRAD
  (p3 ph3):f2
  d0
  (p2 ph2)
  d0
  p16:gp1*EA
  d16
  (p24:sp7 ph4):f2
  DELTA
  p16:gp1*-1*EA
  d16 pl2:f2
  (p3 ph4):f2
  4u
  p16:gp2
  DELTA1 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 
     F1EA(calgrad(EA), caldel(d0, +in0) & calph(ph3, +180) & calph(ph5, +180) & calph(ph31, +180))
exit 
  

ph1=0
ph2=0 0 2 2
ph3=0 2
ph4=0 0 0 0 2 2 2 2
ph5=0
ph31=0 2 0 2 2 0 2 0


;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;sp3: f2 channel - shaped pulse 180 degree
;spnam3: Crp60,0.5,20.1
;sp7: f2 channel - shaped pulse 180 degree for refocussing
;spnam7: Crp60comp.4
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p14: f2 channel - 180 degree shaped pulse for inversion
;     = 500usec for Crp60,0.5,20.1
;p16: homospoil/gradient pulse                         [1 msec]
;p24: f2 channel - 180 degree shaped pulse for refocussing
;     = 2msec for Crp60comp.4
;d0 : incremented delay (2D)                           [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J(XH))
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(XH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;ns: 2 * n
;ds: 16
;td1: number of experiments
;FnMODE: echo-antiecho
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    gp 1 : gp 2
;                         80 : 40.2

;for z-only gradients:
;gpz1: 80%
;gpz2: 40.2%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100



;$Id: hmqcetgp.2,v 1.8 2012/01/31 17:49:23 ber Exp $
