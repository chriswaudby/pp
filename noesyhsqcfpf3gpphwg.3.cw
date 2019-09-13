;NOESY-15N HSQC (derived from hsqcfpf3gpphwg.3.2.jk)
; incorporating flip-back during mixing time using radiation damping
; set d21 = radiation damping period (10-50 ms)
; 45/180 phase correction in t1
;For 15N-labelled samples only (no 13C decoupling)
;With water flipback
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>

  
#include <Avance.incl>
#include <Delay.incl>
#include <Grad.incl>


"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"
;"d13=4u"
"d26=1s/(cnst4*4)"


"in0=inf1/2"
"in10=inf2/4"

"DELTA1=d26-p16-d16-4u"
"DELTA2=d26-p19-d16-p11-d12-4u"
"DELTA3=d26-p19-d16-p11-d12-8u-de+0.63662*p1"

"TAU=d8-p16-d16-p11-20u"


"d0=in0/2-0.6366*p1"
# ifdef N_SINGLEDWELL
    "d10=in10-0.5*(10u+p1+0.63662*p21)"
# else
    "d10=in10/2-0.5*(10u+p1+0.63662*p21)"
# endif /*N_SINGLEDWELL*/

"acqt0=0"

aqseq 321


1 ze 
  d11 pl16:f3
2 d11 do:f3

  d1

3 d12 pl1:f1 pl3:f3
  50u UNBLKGRAD

  ;(p21 ph1):f3
  ;4u
  ;p16:gp0
  ;d16
  ;(p21 ph2):f3
  ;4u
  ;p16:gp0*0.7
  ;d16*2

# ifndef NH_PLANE
  ; 1H evolution
  10u pl1:f1

  ; half-dwell = 1.2*p1+p22+8u+2*d0
  ; => d0 = in1/2 - 0.6366*p1 - p21 - 4u
  "d20=d0-p21-4u"
  if "d20 < 4u"
  {
    (p1 ph10):f1
    d0
    d0
    (p1 ph11):f1
  }
  else
  {
    (p1 ph10):f1
    d20 gron1
    4u groff
    (p22 ph5):f3
    d20 gron1*-1
    4u groff
    (p1 ph11):f1
  }
  ; NOESY mixing period
  "TAU=d8*0.5-d21*0.5-p16*1.5-d16*0.5-4u"
  p16:gp2
  TAU gron3
  4u groff
  TAU gron3*-1
  4u groff
  p16:gp2*-1

  d21 ; radiation damping time [10 ms]

  p16:gp4
  d16
;goto 999

# endif /*NH_PLANE*/
 

  (p1 ph1)
  4u
  p16:gp5
  d16
  DELTA1
  (center (p2 ph2) (p22 ph6):f3 )
  DELTA1
  4u
  p16:gp5
  d16
  (p1 ph2) 

  4u pl0:f1
  (p11:sp1 ph8:r):f1	; flipdown(x), -y -> -z
  4u
  p16:gp6
  d16 pl1:f1

# ifdef HH_PLANE

    (p21 ph3):f3
    2u
    (p21 ph4):f3
    2u
    (p2 ph5):f1

# else

    (p21 ph3):f3
    2u
    d10 gron9
    d10 gron9*-1
    8u groff
    (p2 ph5):f1
    2u
    d10 gron9
    d10 gron9*-1
    8u groff
    (p21 ph4):f3

# endif /*HH_PLANE*/

  4u
  p16:gp8
  d16 pl0:f1
  (p11:sp1 ph7:r):f1	; flipdown(-x), z -> y
  4u
  4u pl1:f1

  (p1 ph1)
  4u
  p19:gp7
  d16 pl0:f1
  DELTA2
  (p11:sp1 ph7:r):f1	; flipdown(-x), z -> y
  d12 pl1:f1
  (center (p2 ph1) (p22 ph1):f3 )
  d12 pl0:f1
  (p11:sp11 ph9:r):f1	; flipback(-x), -y -> z
  4u
  p19:gp7
999  d16 pl16:f3
  DELTA3
  4u BLKGRAD

  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
	F1PH(rp3 & rp6 & rd10 & ip11 & ip21, id0)
	F2PH(ip3 & ip6, id10)
exit 
  

ph1=0
ph2=1
ph6=0
ph7=2
ph8=0
ph9=2

# ifndef NH_PLANE

ph3=0 2
ph4=0 0 0 0 2 2 2 2
ph5=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph10=(8) 1 1 5 5
ph11=0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
ph20=2 2 0 0
ph21=2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
ph31=0 2 2 0 2 0 0 2 0 2 2 0 2 0 0 2
     2 0 0 2 0 2 2 0 2 0 0 2 0 2 2 0

# else

ph3=0 2
ph4=0 0 2 2
ph5=0 0 0 0 2 2 2 2
ph31=0 2 2 0

# endif /*NH_PLANE*/


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl16: f3 channel - power level for CPD/BB decoupling
;sp1: f1 channel - shaped pulse  90 degree (flip-down)
;sp11: f1 channel - shaped pulse  90 degree (flip-back)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse
;p19: second homospoil/gradient pulse
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;d10 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d8 : mixing time
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21: radiation damping delay [10 msec]
;d26 : 1/(4J)NH
;cnst4: = J(NH)
;inf1: 1/SW(H) = 4 * DW(H)
;inf2: 1/SW(N) = 4 * DW(N)
;in0: 1/(2 * SW(H)) = DW(H)
;in10: 1/(2 * SW(N)) = DW(N)
;nd0: 4
;nd10: 4
;NS: 4 * n
;DS: 16
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI, TPPI, States or QSEQ in F1
;FnMODE: States-TPPI, TPPI, States or QSEQ in F2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz9: 1-2%
;gpz0: 47%
;gpz1: 2%
;gpz2: 29%
;gpz3: 2%
;gpz4: 35%
;gpz5: 7%
;gpz6: 17%
;gpz7: 53%
;gpz8: 13%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100

