;noesygpph.cw
;with addition of a purge pulse before d1
;avance-version (12/01/11)
;2D homonuclear correlation via dipolar coupling 
;dipolar coupling may be due to noe or chemical exchange.
;phase sensitive
;with gradient pulses in mixing time
;
;J. Jeener, B.H. Meier, P. Bachmann & R.R. Ernst, J. Chem. Phys. 71,
;   4546-4553 (1979)
;R. Wagner & S. Berger, J. Magn. Reson. 123 A, 119-121 (1996)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<delay> tmix=<$VDLIST>
"l1=0"

"p2=p1*2"
"d12=1ms"

"in0=inf2"

"d0=in0/2-p1*4/3.1416"


;"TAU=d8*0.5-p16-d16-50u"


"acqt0=-p1*2/3.1416"

aqseq 312

1 ze
2 d12
  "TAU=tmix[l1]*0.5-p16-d16-50u"

  20u pl11:f1
  (2mp ph10):f1
  20u
  (3mp ph11):f1
  20u pl1:f1

  d1
3 p1 ph1
  d0
  p1 ph2
  TAU 
  50u UNBLKGRAD
  p16:gp1
  d16
  3u
  (p2 ph4):f1
  3u
  p16:gp1*-1
  d16
  50u BLKGRAD
  TAU 
  p1 ph3
  go=2 ph31
  d12 mc #0 to 2
    F1QF(calclc(l1,1))
    F2PH(calph(ph1, +90), caldel(d0, +in0))
exit


ph1=0 2 
ph2=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph3=0 0 2 2 1 1 3 3
ph4=0
ph10=0
ph11=1
ph31=0 2 2 0 1 3 3 1 2 0 0 2 3 1 1 3


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse                       [1 msec]
;d0 : incremented delay (2D)
;d1 : relaxation delay; 1-5 * T1
;d8 : mixing time
;d16: delay for homospoil/gradient recovery
;inf1: 1/SW = 2 * DW
;in0: 1/(1 * SW) = 2 * DW
;nd0: 1
;ns: 2 * n
;ds: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ


;use gradient ratio:    gp 1
;                         40

;for z-only gradients:
;gpz1: 40%

;use gradient files:   
;gpnam1: SMSQ10.100


;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1



;$Id: noesygpph,v 1.11 2012/01/31 17:49:27 ber Exp $
