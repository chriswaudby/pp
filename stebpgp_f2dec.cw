;with 13C CPD during acquisition
;With added presaturation (during recycle delay only)
;Zero first-order phase correction
;
;stebpgp1s
;avance-version (07/05/08)
;2D sequence for diffusion measurement using stimulated echo
;using bipolar gradient pulses for diffusion
;using 1 spoil gradient
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple_d>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


define list<gradient> diff=<Difframp>


"p2=p1*2"

"d11=30m"
"d12=20u"
"d13=4u"


"TAU=8u+de-0.6366*p1"
"DELTA1=d20-p1*2-p2-p30*2-d16*3-p19-TAU"

"acqt0=de"

1 ze
  4u pl12:f2
  d12 pl1:f1
2 d1 
  50u UNBLKGRAD
  p1 ph1
  p30:gp6*diff
  d16
  p2 ph2
  p30:gp6*-1*diff
  d16
  p1 ph3
  p19:gp7
  d16 
  DELTA1
  p1 ph4
  TAU
  p30:gp6*diff
  d16
  p2 ph2
  p30:gp6*-1*diff
  d16 pl12:f2
  8u BLKGRAD
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 F1QF(igrad diff)
exit


ph1= 0
ph2= 0 0 0 0  2 2 2 2
ph3= 0 0 0 0  0 0 0 0  2 2 2 2  2 2 2 2
ph4= 0 1 2 3
ph29=0
ph31=0 3 2 1  0 3 2 1  2 1 0 3  2 1 0 3


;pl1: f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p19: gradient pulse 2 (spoil gradient)
;p30: gradient pulse (little DELTA * 0.5)
;d1 : relaxation delay; 1-5 * T1
;d16: delay for gradient recovery
;d20: diffusion time (big DELTA)
;NS : 8 * n
;DS : 4 * m
;td1: number of experiments
;FnMODE: QF
;        use xf2 and DOSY processing


;use gradient ratio:    gp 6 : gp 7
;                       100  : -17.13

;for z-only gradients:
;gpz6: 100%
;gpz7: -17.13% (spoil)

;use gradient files:   
;gpnam6: SINE.100
;gpnam7: SINE.100

;use AU-program dosy to calculate gradient-file Difframp



;$Id: stebpgp1s,v 1.4.6.1 2007/05/09 09:36:59 ber Exp $
