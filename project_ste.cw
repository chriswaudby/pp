; project_ste.cw
; PROJECT-STE T2-D correlation sequence
;
; based on J.A. Aguilar, M. Nilsson, G. Bodenhausen and G. Morris
; Chem. Commun. (2012) 48 811-813
;
; run as pseudo-3D
; with presaturation
;
; Chris Waudby
; July 2018

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

prosol relations=<triple>

;define list<gradient> diff=<Difframp>
define list<gradient> diff={0.089 0.096 0.121 0.164 0.214 0.271 0.332 0.396 0.463 0.532 0.604 0.675 0.750 0.827 0.904 0.982}
;define list<loopcounter> ncyc = {0 1 2 4 6 8 10 13 16 19 22 26 29 33 37 41}
define list<loopcounter> ncyc = {0 41 1 37 2 33 4 29 6 26 8 22 10 19 13 16}

"p2=p1*2"
"d2=250u"
"d3=4*d2"  ; T2 time per loop
"d12=20u"
"DELTA1=d2-p1"
"DELTA2=d2-1.5*p1"
"TAU=d20-p30*2-p19-d16*3-p11-p1*4-8u"
"acqt0=-0.6366*p1"

; acquisition order = D then T2
aqseq 321

1 ze
2 30m pl9:f1
  d1 cw:f1 ph29
  4u do:f1
  d12 pl1:f1 UNBLKGRAD

  ; start PROJECT echos
  p1 ph11

3 DELTA1
  p2 ph2
  DELTA2
  p1 ph12
  DELTA2
  p2 ph2
  DELTA1
  lo to 3 times ncyc

  ; STE encoding
  p30:gp6*diff
  d16
  p2 ph14
  p30:gp6*-1*diff
  d16
  p1 ph15

  ; diffusion delay
  4u
  p19:gp7
  d16 ;pl9:f1
  TAU ;cw:f1 ph1
  4u ;do:f1

  ; water flip-down and STE decoding
  (p11:sp1 ph23:r)
  4u pl1:f1
  p1 ph13
  p30:gp6*diff
  d16
  p2 ph1
  p30:gp6*-1*diff
  d16
  4u BLKGRAD

  ; acquisition
  go=2 ph31
  30m mc #0 to 2
      F2QF(igrad diff)
      F1QF(ncyc.inc)
exit

ph1 =0
ph2 =1
ph11=0 2
ph12=1 1 3 3
ph13=0 0 0 0  3 3 3 3  2 2 2 2  1 1 1 1
ph23=2 2 2 2  1 1 1 1  0 0 0 0  3 3 3 3
ph14=0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0
     2 2 2 2  2 2 2 2  2 2 2 2  2 2 2 2
ph15=0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  
     0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  
     2 2 2 2  2 2 2 2  2 2 2 2  2 2 2 2  
     2 2 2 2  2 2 2 2  2 2 2 2  2 2 2 2  
ph29=0
ph31=0 2 0 2  1 3 1 3  2 0 2 0  3 1 3 1
     0 2 0 2  1 3 1 3  2 0 2 0  3 1 3 1
     2 0 2 0  3 1 3 1  0 2 0 2  1 3 1 3
     2 0 2 0  3 1 3 1  0 2 0 2  1 3 1 3

;pl1: f1 channel - power level for pulse (default)
;pl9: f1 channel - power level for presaturation
;p1 : f1 channel - 90 degree excitation pulse
;p2 : f1 channel - 180 degree refocusing pulse
;p19: gradient pulse 2 (spoil gradient)
;p30: gradient pulse (little DELTA * 0.5) [2.75ms]
;d1 : relaxation delay; 1-5 * T1
;d2 : tau (T2 time = 4*d2*vc)
;d3 : total T2 time per loop
;d12: delay for power switching [20 usec]
;d16: delay for gradient recovery
;d20: diffusion time (big DELTA)
;ds : 16 
;ns : 32*n

;for z-only gradients:
;gpz6: 100%
;gpz7: -17.13% (spoil)

;use gradient files:   
;gpnam6: SMSQ10.100
;gpnam7: SINE.100

;use AU-program dosy to calculate gradient-file Difframp

