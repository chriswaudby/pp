;cosyetgp
; with continuous cpd in f2
; corrected for zero first order phase correction
; using baseopt
;avance-version (12/01/11)
;2D homonuclear shift correlation
;phase sensitive using Echo/Antiecho-TPPI gradient selection
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
"d11=30m"


"in0=inf1"

"d0=3u"


"DELTA=p16+d16+d0"
"DELTA1=p16+d16+8u"
"acqt0=0"


1 ze
  d11 pl12:f2
  d11 
2 d11 do:f2
  d1
  50u UNBLKGRAD
  50u cpd2:f2

  p1 ph1

  DELTA

  p2 ph2

  d0
  p16:gp1*EA
  d16

  p1 ph3

  DELTA1

  p2 ph2

  4u
  p16:gp2
  d16
  4u BLKGRAD

  go=2 ph31
  d11 do:f2 mc #0 to 2
      F1EA(calgrad(EA), caldel(d0, +in0) & calph(ph1, +180) & calph(ph31, +180))
  d11 do:f2
exit


ph1=0 2 2 0 1 3 3 1
ph2=0 2 0 2 1 3 1 3
ph3=0 2 0 2 1 3 1 3
ph31=0 2 2 0 1 3 3 1


;pl1 : f1 channel - power level for pulse (default)
;p0 : f1 channel -  20 to 90 degree high power pulse
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d16: delay for homospoil/gradient recovery
;inf1: 1/SW = 2 * DW
;in0: 1/(1 * SW) = 2 * DW
;nd0: 1
;ns: 4 * n
;ds: 16
;td1: number of experiments
;FnMODE: echo-antiecho


;use gradient ratio:	gp 1 : gp 2
; 		          30 :   30

;for z-only gradients:
;gpz1: 30%
;gpz2: 30%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100



