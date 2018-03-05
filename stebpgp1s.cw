;stebpgp1s
;using baseopt
;avance-version (09/04/17)
;2D sequence for diffusion measurement using stimulated echo
;using bipolar gradient pulses for diffusion
;using 1 spoil gradient
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


define list<gradient> diff=<Difframp>


"p2=p1*2"


"DELTA1=d20-p1*2-p2-p30*2-d16*3-p19"
"DELTA2=d16-4u+p1*0.6366"
"acqt0=0"
baseopt_echo

1 ze
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
  p30:gp6*diff
  d16
  p2 ph2
  p30:gp6*-1*diff
  DELTA2
  4u BLKGRAD

  go=2 ph31 
  d1 mc #0 to 2 F1QF(calgrad(diff))

exit


ph1= 0
ph2= 0 0 0 0  2 2 2 2
ph3= 0 0 0 0  0 0 0 0  2 2 2 2  2 2 2 2
ph4= 0 1 2 3
ph31=0 3 2 1  0 3 2 1  2 1 0 3  2 1 0 3


;pl1: f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p19: gradient pulse 2 (spoil gradient)
;p30: gradient pulse (little DELTA * 0.5)
;d1 : relaxation delay; 1-5 * T1
;d16: delay for gradient recovery
;d20: diffusion time (big DELTA)
;ns : 8 * n
;ds : 4 * m
;td1: number of experiments
;FnMODE: QF
;        use xf2 and DOSY processing


;use gradient ratio:    gp 6 : gp 7
;                       100  : -17.13

;for z-only gradients:
;gpz6: 100%
;gpz7: -17.13% (spoil)

;use gradient files:   
;gpnam6: SMSQ10.100
;gpnam7: SMSQ10.100

;use AU-program dosy to calculate gradient-file Difframp



;$Id: stebpgp1s,v 1.7.2.1 2013/08/30 09:43:35 ber Exp $
