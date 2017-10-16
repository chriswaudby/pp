;Water-LOGSY sequence for ligand-screening
;Water inversion using e-PHOGSY element with water-selective ReBURP pulse
;Suppression of auto-relaxed signals via phase-cycling of ReBURP pulse
;  and 180deg pulse during mixing time
;Gradients during mixing time to prevent radiation damping
;Water returned to z at end of sequence
;Water suppression using excitation sculpting with gradients
;Full EXORCYCLE on soft 180--hard 180 pulse pairs in water-DPFGSE element
;With spin-lock period for suppression of receptor signals
;John K, Jan 2012

;C. Dalvit & U. Hommel, J. Magn. Reson. Ser. B 109, 334-338 (1995)
;C. Dalvit et al, J. Biomol. NMR 18, 65-68 (2000)
;C. Dalvit et al, J. Biomol. NMR 21, 349-359 (2001)

;T.-L. Hwang & A.J. Shaka, J. Magn. Reson. Ser. A 112 275-279 (1995)
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<logsy>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"d11=30m"
"d12=20u"
"d13=4u"

"d18=0"		; comment out to allow manual setting of d18

"TAU=0.5*(d8-p2-p16-d16-100u)"

"TAU1=de+0.63662*p1+50u"

"l0=1"

"p29=d29"

1 ze
2 d11
  d12 pl1:f1 BLKGRAD
  d1
  50u UNBLKGRAD

; e-PHOGSY

  (p1 ph1):f1
  p19:gp1
  d16 pl0:f1
  d18
  (p32:sp22 ph2:r):f1
  d18
  p19:gp1
  d16 pl1:f1
  (p1 ph3):f1

; Mixing time

  p16:gp2
  d16
  TAU gron0
  50u groff
  (p2 ph1):f1
  10u
  TAU gron0*-1
  40u groff

; Spin-lock

  (p1 ph10):f1
  4u pl10:f1
  (p29 ph11):f1
  4u pl1:f1
  (p1 ph12):f1
  p16:gp5
  d16

; Water flipback/down

  d12 pl0:f1
  if "l0%2 == 1"
     {
     (p11:sp1 ph4:r):f1		; flipback (phase +x), -z -> y
     }
  else
     {
     (p11:sp21 ph14:r):f1	; flipdown (phase -x), +z -> y
     }
  4u iu0
  d12 pl1:f1

; water-DPFGSE

  (p1 ph5):f1
  
  50u
  p16:gp3
  d16 pl0:f1
  (p12:sp11 ph6:r):f1		; 180deg flipdown
  d12 pl1:f1
  (p2 ph7):f1
  p16:gp3
  d16 

  TAU1

  p16:gp4
  d16 pl0:f1
  (p12:sp11 ph8:r):f1		; 180deg flipdown
  d12 pl1:f1
  (p2 ph9):f1
  p16:gp4
  d16

  go=2 ph31
  d11 mc #0 to 2 F0(zd)
  d13 BLKGRAD

exit


ph1= 0
ph2= 0 1 2 3
ph3= 0
ph4= 0 2 0 2	; 90 deg flipback (0 _ 0 _)
ph14=0 2 0 2	; 90 deg flipdown (_ 2 _ 2)
ph5= 0

ph6= 0 0 0 0  2 2 2 2  0 0 0 0  2 2 2 2
     1 1 1 1  3 3 3 3  1 1 1 1  3 3 3 3
ph7= 2 2 2 2  0 0 0 0  2 2 2 2  0 0 0 0
     3 3 3 3  1 1 1 1  3 3 3 3  1 1 1 1

ph8= 0 0 0 0  0 0 0 0  2 2 2 2  2 2 2 2 
     0 0 0 0  0 0 0 0  2 2 2 2  2 2 2 2
     1 1 1 1  1 1 1 1  3 3 3 3  3 3 3 3
     1 1 1 1  1 1 1 1  3 3 3 3  3 3 3 3
ph9= 2 2 2 2  2 2 2 2  0 0 0 0  0 0 0 0 
     2 2 2 2  2 2 2 2  0 0 0 0  0 0 0 0
     3 3 3 3  3 3 3 3  1 1 1 1  1 1 1 1
     3 3 3 3  3 3 3 3  1 1 1 1  1 1 1 1

ph10=1
ph11=0
ph12=3

ph31=0 2 0 2  0 2 0 2  0 2 0 2  0 2 0 2
     2 0 2 0  2 0 2 0  2 0 2 0  2 0 2 0
     2 0 2 0  2 0 2 0  2 0 2 0  2 0 2 0
     0 2 0 2  0 2 0 2  0 2 0 2  0 2 0 2
     
     


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl10 : f1 channel - power level for spin-lock
;sp1 : f1 channel - shaped pulse 90 degree [Eburp flipback]
;sp11 : f1 channel - shaped pulse 180 degree [Eburp flipdown]
;sp21 : f1 channel - shaped pulse 90 degree [flipdown]
;sp22 : f1 channel - shaped pulse 180 degree [e-PHOGSY]
;spnam1 : Eburp.1000
;spnam11 : Squa100.1000
;spnam21 : Eburp.1000
;spnam22 : Reburp.1000
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p11: f1 channel -  90 degree shaped pulse (Eburp.1000)      [4 msec]
;p12: f1 channel - 180 degree shaped pulse (Squa100.1000)   [2 msec]
;p32: f1 channel - 180 degree shaped pulse (Reburp.1000)    [5-25 msec]
;p16: homospoil/gradient pulse			 [1 msec]
;p19: gradient pulse in e-PHOGSY		 [1 msec]
;d1 : relaxation delay; 1-5 * T1		 [2-5 sec]
;d8 : mixing time				 [1-2 sec]
;d11: delay for disk I/O                         [30 msec]
;d12: delay for power switching                  [20 usec]
;d13: short delay                                [4 usec]
;d16: delay for homospoil/gradient recovery      [200 usec]
;d18: delay for attenuation of receptor resonances at water shift
;d29: spin-lock time				 [20-200 msec]
;NS: 16 * n, total number of scans: NS * TD0
;DS: 8


;for z-only gradients:
;gpz0: 2%
;gpz1: 17%
;gpz2: 43%
;gpz3: 31%
;gpz4: 13%
;gpz5: 37%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100

