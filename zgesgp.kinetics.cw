;zgesgp.kinetics.cw
;avance-version (07/10/04)
;pseudo-2D sequence
;water suppression using excitation sculpting with gradients
;T.-L. Hwang & A.J. Shaka, J. Magn. Reson.,
;   Series A 112 275-279 (1995)
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"d12=20u"


"TAU=de+p1*2/3.1416+50u"
"acqt0=de"

1 ze
2 30m
3 4u
#ifdef PRESAT
  d12 pl9:f1 BLKGRAD
  d1 cw:f1
  4u do:f1
  4u pl1:f1
#else
  d12 pl1:f1 BLKGRAD
  d1
  8u
#endif
  p1 ph1
  
  50u UNBLKGRAD
  p16:gp1
  d16 pl0:f1
  (p12:sp11 ph2:r):f1
  4u
  d12 pl1:f1

  p2 ph3

  4u
  p16:gp1
  d16 
  TAU
  p16:gp2
  d16 pl0:f1
  (p12:sp11 ph4:r):f1
  4u
  d12 pl1:f1

  p2 ph5

  4u
  p16:gp2
  d16

  go=2 ph31
  30m wr #0 if #0 zd
  lo to 3 times td1

  4u BLKGRAD


exit


ph1=0
ph2= 0 1 0 1  2 3 2 3  0 1 0 1  2 3 2 3
ph3= 2 3 2 3  0 1 0 1  2 3 2 3  0 1 0 1 
ph4= 0 0 1 1  0 0 1 1  2 2 3 3  2 2 3 3 
ph5= 2 2 3 3  2 2 3 3  0 0 1 1  0 0 1 1 
ph31=0 2 2 0  0 2 2 0  0 2 2 0  0 2 2 0


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;sp1 : f1 channel - shaped pulse 180 degree
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p12: f1 channel - 180 degree shaped pulse (Squa100.1000)   [2 msec]
;p16: homospoil/gradient pulse
;d1 : relaxation delay; 1-5 * T1
;d12: delay for power switching                             [20 usec]
;d16: delay for homospoil/gradient recovery
;NS: 8 * n, total number of scans: NS * TD0
;DS: 4


;use gradient ratio:    gp 1 : gp 2
;                         31 :   11

;for z-only gradients:
;gpz1: 31%
;gpz2: 11%

;use gradient files:
;gpnam1: SINE.100
;gpnam2: SINE.100



;$Id: zgesgp,v 1.5.6.1 2007/10/04 16:52:07 ber Exp $
