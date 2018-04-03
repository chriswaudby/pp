;pseudo-3D for kinetic measurements
; Delays in final spin-echo adjusted to give zero first-order phase corr.
; Use for flat baselines (with baseopt, Chris Waudby October 2016)
;
;stebpgp1s19
;avance-version (07/05/08)
;2D sequence for diffusion measurement using stimulated echo
;using bipolar gradient pulses for diffusion
;using 1 spoil gradient
;water suppression using 3-9-19 pulse sequence with gradients
;with added presat water suppression
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


;define list<gradient> diff=<Difframp>
;define list<gradient> diff={0.05 0.179 0.307 0.436 0.564 0.693 0.821 0.95}
define list<gradient> diff={0.15 0.95}

"p2=p1*2"
"d11=30m"

"DELTA1=d20-p1*2-p2-p30*2-d16*3-p19"

"TAU=0.6366*p1+8u"

"acqt0=0"
aqseq 321

1 ze
2 d11
  4u pl9:f1
  d1 cw:f1 ph1
  4u do:f1

  50u pl1:f1 UNBLKGRAD
  p1 ph1
  p30:gp6*diff
  d16
  p2 ph2
  p30:gp6*-1*diff
  d16
  p1 ph3
  p19:gp7

  d16 pl9:f1
  DELTA1 cw:f1 ph1
  1u do:f1
  4u pl1:f1
;  d16  
;  DELTA1
  p1 ph4
  p30:gp6*diff
  d16
  p2 ph2
  p30:gp6*-1*diff
  d16
  p16:gp1
  d16 pl18:f1
  TAU
  p27*0.231 ph5
  d19*2
  p27*0.692 ph5
  d19*2
  p27*1.462 ph5
  d19*2
  p27*1.462 ph6
  d19*2
  p27*0.692 ph6
  d19*2
  p0*0.231 ph6
  4u
  p16:gp1
  d16
  4u BLKGRAD
  go=2 ph31 
  d11 mc #0 to 2
    F2QF(igrad diff)
    F1QF(rgrad diff)
exit


ph1= 0
ph2= 0 0 0 0  2 2 2 2
ph3= 0 0 0 0  0 0 0 0  2 2 2 2  2 2 2 2
ph4= 0 1 2 3
ph5= 0
ph6= 2
ph31=0 1 2 3  0 1 2 3  2 3 0 1  2 3 0 1


;pl1: f1 channel - power level for pulse (default)
;pl9: f1 channel - power level for presaturation
;pl18: f1 channel - power level for 3-9-19-pulse (watergate)
;p0 : f1 channel -  90 degree pulse at pl18
;                      use for fine adjustment
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p16: gradient pulse (WATERGATE)
;p19: gradient pulse 2 (spoil gradient)
;p27: f1 channel -  90 degree pulse at pl18
;p30: gradient pulse (little DELTA * 0.5)
;d1 : relaxation delay; 1-5 * T1
;d16: delay for gradient recovery
;d19: delay for binomial water suppression
;     d19 = (1/(2*d)), d = distance of next null (in Hz) [106 us at 600 MHz]
;d20: diffusion time (big DELTA)
;NS : 8 * n
;DS : 4 * m
;td1: number of experiments
;FnMODE: QF
;        use xf2 and DOSY processing


;use gradient ratio:    gp 1 : gp 6 : gp7
;                        -20 :  100 : -17.13

;for z-only gradients:
;gpz1: -20%
;gpz6: 100%
;gpz7: -17.13% (spoil)

;use gradient files:   
;gpnam1: SINE.100
;gpnam6: SINE.100
;gpnam7: SINE.100

;use AU-program dosy to calculate gradient-file Difframp



;$Id: stebpgp1s19,v 1.4.6.1 2007/05/09 09:36:59 ber Exp $
