;TRACT (TROSY for rotational correlation times)
;Lee, Hilty, Wider & Wuethrich, J. Magn. Reson. 178, 72-76 (2006)
;Measures decay of TROSY line by default
;Use ZGOPTN -DANTITROSY to measure decay of anti-TROSY line

;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=

prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p22=p21*2"

"d11=30m"
"d12=20u"
;"d13=4u"
"d26=1s/(cnst4*4)"


"DELTA1=d26-p19-d16-p11-2*d12"
"DELTA2=d26-p19-d16-p11-2*d12-p21"


1 ze 
2 d11
3 d11 pl1:f1 pl3:f3
  d12
    "TAU=vd/2-p21"
  d1
  50u UNBLKGRAD

  (p21 ph1):f3
  4u
  p16:gp1
  d16*2

  (p1 ph1):f1

  d12 pl0:f1
  (p11:sp11 ph3:r):f1	; flipback(-x): -y -> +z
  d12 pl1:f1
  DELTA1
  p19:gp2
  d16

  ( center (p2 ph2):f1 (p22 ph1):f3 )

  DELTA1
  p19:gp2
  d16
  d12 pl0:f1
  (p11:sp11 ph2:r):f1	; flipback(+y): -z -> -x
  d12 pl1:f1

  ( center (p1 ph2):f1 (p21 ph11):f3 )

  TAU
  (p22 ph1):f3
  TAU

  (p1 ph4):f1

  d12 pl0:f1
  (p11:sp11 ph2:r):f1	; flipback(+y): -x -> +z
  d12 pl1:f1
  DELTA1
  p19:gp3
  d16

  ( center (p2 ph1):f1 (p22 ph4):f3 )

  DELTA1
  p19:gp3
  d16
  d12 pl0:f1
  (p11:sp11 ph1:r):f1	; flipback(+x): -z -> +y
  d12 pl1:f1

  ( center (p1 ph1):f1 (p21 ph4):f3 )

  DELTA1
  p19:gp4
  d16
  d12 pl0:f1
  (p11:sp1 ph13:r):f1	; flipdown(-x): +z -> +y
  d12 pl1:f1

  ( center (p2 ph1):f1 (p22 ph1):f3 )

  d12 pl0:f1
  (p11:sp11 ph3:r):f1	; flipback(-x): -y -> +z
  d12 pl1:f1
  p19:gp4
  d16 BLKGRAD
  DELTA2

  (p21 ph12):f3

  go=2 ph31
  d11 mc #0 to 2 F1QF(ivd)
exit 
  

ph1= 0
ph2= 1
ph3= 2
ph13=2		; flipdown only
ph4= 3

# ifdef ANTITROSY
ph11=1 3 2 0
ph12=2
# else 
ph11=1 3 0 2
ph12=0
# endif /*ANTITROSY*/

ph31=1 3 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;sp1 : f1 channel - shaped pulse  90 degree [flipdown]
;sp11: f1 channel - shaped pulse  90 degree [flipback]
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p11: f1 channel -  90 degree shaped pulse [water selective]
;p16: homospoil/gradient pulse                       [1 msec]
;p19: homospoil/gradient pulse                       [500 usec] 
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power-switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d26 : 1/(4J)YH
;cnst4: = J(YH)
;NS: 4 * n
;DS: 4
;td1: number of experiments
;FnMODE: QF


;for z-only gradients:
;gpz1: 80%
;gpz2: 17%
;gpz3: 13%
;gpz4: 53%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100

