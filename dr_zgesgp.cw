;dr_zgesgp.cw
;avance-version (17/02/02)
;1D sequence
;for dual receive on f1 and f2
; 19F 1D (parent) on f1
; 1H excitation sculpting (child) on f2
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <De.incl>
#include <Grad.incl>
#include <Delay.incl>

"p4=p3*2" ; 1H hard 180

"TAU=50u+p3*0.6366+p1*0.3634"
"acqt0=-p1*2/3.1416"  ; optimise baseopt for 19F


1 ze1
  ze2
  30m
2 d12 BLKGRAD
  d1
  d12 pl1:f1
  d12 pl2:f2

  (p3 ph1):f2
  50u UNBLKGRAD
  p16:gp1
  d16
  (p12:sp1 ph2:r):f2
  4u
  d12 pl2:f2

  (p4 ph3):f2

  4u
  p16:gp1
  d16
  TAU
  p16:gp2
  d16
  (p12:sp1 ph4:r):f2
  4u
  d12 pl2:f2

  (p4 ph5):f2

  4u
  p16:gp2
  d16

  ; 19F excitation
  (p1 ph10):f1

  ; start dual receiver acquisition
  ACQ_START1(ph30,ph11) ACQ_START2(ph30,ph31)
  0.1u DWELL_GEN1 DWELL_GEN2
  (aq1) (aq2)
  eoscnp2
  rcyc=2 
  30m wr1 #0
  30m wr2 #1
exit

; reference phase
ph30=0
; 19F phases
ph10=0 2 2 0 1 3 3 1
ph11=0 2 2 0 1 3 3 1
; 1H phases
ph1=0
ph2= 0 1 0 1  2 3 2 3  0 1 0 1  2 3 2 3
ph3= 2 3 2 3  0 1 0 1  2 3 2 3  0 1 0 1
ph4= 0 0 1 1  0 0 1 1  2 2 3 3  2 2 3 3
ph5= 2 2 3 3  2 2 3 3  0 0 1 1  0 0 1 1
ph31=0 2 2 0  0 2 2 0  0 2 2 0  0 2 2 0


;pl2 : f2 channel (1H) - power level for pulse (default)
;pl1 : f1 channel (19F) - power level for pulse (default)
;sp1 : f2 channel - shaped pulse 180 degree
;p1 : f1 channel (19F) -  90 degree high power pulse
;p3 : f2 channel (1H) -  90 degree high power pulse
;p4 : f2 channel (1H) - 180 degree high power pulse
;p12: f2 channel - 180 degree shaped pulse (Squa100.1000)   [2 msec]
;p16: homospoil/gradient pulse
;d1 : relaxation delay; 1-5 * T1
;d12: delay for power switching                             [20 usec]
;d16: delay for homospoil/gradient recovery
;ns: 8 * n

;use gradient ratio:    gp 1 : gp 2
;                         31 :   11

;for z-only gradients:
;gpz1: 31%
;gpz2: 11%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100

