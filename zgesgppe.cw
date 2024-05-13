;zgesgppe.cw
;avance-version (13/08/01)
;1D sequence
;water suppression using excitation sculpting with gradients
;   using perfect echo
;
;(R.W. Adams, C.M. Holroyd, J.A. Aguilar, M. Nilsson & G.A. Morris,
;   Chem. Commun. 49, 358-360 (2013))
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

/*  water suppression (p12, sp11)   */
"p12=2000u"
"spw1=plw1*(pow(p1*2/p12,2))" /* power level for Squa100.1000  */
"spoffs1=0"
"spoal1=0.5"


"DELTA1=p12+p16+d16+p2/2+de/2+p1/PI+12u"
"TAU=de+p1*2/PI"


"acqt0=0"
baseopt_echo


1 ze
2 30m
  d12 pl1:f1 BLKGRAD
  d1
  50u UNBLKGRAD

  (p1 ph1)

  p16:gp3
  d16
  DELTA1
  (p2 ph7)
  DELTA1
  p16:gp3
  d16

  (p1 ph6)
  
  p16:gp1
  d16
  (p12:sp1 ph2:r):f1
  4u
  4u pl1:f1

  p2 ph3

  4u
  p16:gp1
  d16 
  TAU
  p16:gp2
  d16
  (p12:sp1 ph4:r):f1
  4u
  4u pl1:f1

  p2 ph5

  4u
  p16:gp2
  d16

  go=2 ph31
  30m mc #0 to 2 F0(zd)
  4u BLKGRAD
exit


ph1=0
ph2=0 1
ph3=2 3
ph4=0 0 1 1
ph5=2 2 3 3
ph6=1
ph7=0
ph31=0 2 2 0 


;pl1 : f1 channel - power level for pulse (default)
;sp1 : f1 channel - shaped pulse 180 degree
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p12: f1 channel - 180 degree shaped pulse (Squa100.1000)   [2 msec]
;p16: homospoil/gradient pulse
;d1 : relaxation delay; 1-5 * T1
;d12: delay for power switching                             [20 usec]
;d16: delay for homospoil/gradient recovery
;ns: 8 * n, total number of scans: NS * TD0
;ds: 4
;spnam11: Squa100.1000

;for z-only gradients:
;gpz1: 31%
;gpz2: 11%
;gpz3: 5%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100



;$Id:$
