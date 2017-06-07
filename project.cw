; project.cw
; PROJECT sequence
;
; J.A. Aguilar, M. Nilsson, G. Bodenhausen and G. Morris
; Chem. Commun. (2012) 48 811-813
;
; run as pseudo-2D
; with presaturation
;
; Chris Waudby
; October 2012

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

prosol relations=<triple>

"p2=p1*2"
"d3=4*d2"  ; T2 time per loop
"d12=20u"
"TAU=50u+de-0.6366*p1"
"DELTA1=d2-p1"
"DELTA2=d2-1.5*p1"

1 ze
2 30m pl9:f1
  d1 cw:f1 ph29
  4u do:f1
  d12 pl1:f1
  p1 ph1
3 DELTA1
  p2 ph2
  DELTA2
  p1 ph3
  DELTA2
  p2 ph2
  DELTA1
  lo to 3 times c
;------------------excitation sculpting   
  p16:gp1
  d16 pl0:f1
  (p12:sp11 ph5:r):f1
  4u
  d12 pl1:f1
  p2 ph6
  4u
  p16:gp1
  d16
  TAU
  p16:gp2
  d16 pl0:f1
  (p12:sp11 ph7:r):f1
  4u
  d12 pl1:f1
  p2 ph8
  4u
  p16:gp2
  d16
  50u BLKGRAD

  go=2 ph31
  30m mc #0 to 2
      F1QF(ivc)
exit

ph1 =0 2
ph2 =1
ph3 =1 1 3 3
ph5 =0 0 0 0 1 1 1 1
ph6 =2 2 2 2 3 3 3 3
ph7 =0 0 0 0 1 1 1 1
ph8 =2 2 2 2 3 3 3 3
ph29=0
ph31=0 2 0 2

;pl1: f1 channel - power level for pulse (default)
;pl9: f1 channel - power level for presaturation
;p1 : f1 channel - 90 degree excitation pulse
;p2 : f1 channel - 180 degree refocusing pulse
;d1 : relaxation delay; 1-5 * T1
;d2 : tau (T2 time = 4*d2*vc)
;d3 : total T2 time per loop
;d12: delay for power switching [20 usec]
;ns : 4*n
