;H(Z/D)QC
;
; simultaneous acquisition of HZQC and HDQC
; separated post-acquisition
; water flip-back
; (90,-180) phase correction
; Crushing equm 15N magnetisation
;
; run as pseudo-3D
; set td1 = 4, QF
; nbl = 4
; ns = 1 (4 step phase cycle stored in F1 dimension)
;
;
; based on Jarvet and Allard (1996) JMR B 112 240-244
;
;K.V. Pervushin, G. Wider, R. Riek & K. Wuethrich, 
;   Proc. Natl. Acad. Sci USA 96, 9607-9612 (1999)
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
"p22=p21*2"
"d11=30m"
"d12=20u"
"d26=1s/(cnst4*4)"


"in0=inf2/2"

"d0=in0/2-p21*0.6366-5u"


"DELTA1=d26-p17-d16"
"DELTA2=d26-p16-d16-p21-p11-4u"
"DELTA3=d26-p16-d16-p21-p11-8u-de"
"acqt0=de"


"l0=1"

"l2=td2/2" ; number of complex points
;"nbl=2"

aqseq 312


1 ze 
2 30u
3 30u
4 d11
5 30u
6 30u

  d1
  d12 pl1:f1 pl3:f3

  ; purge equilibrium Nz
  50u UNBLKGRAD
  (p21 ph1):f3
  p16:gp3
  d16

  ; water flip-down
  if "l0 %2 == 1"
     {
     (p11:sp1 ph2):f1 ; -x
     }
  else
     {
     (p11:sp2 ph1):f1 ; x
     }
  4u pl1:f1

  ; INEPT
  (p1 ph1)
  p17:gp1
  d16
  DELTA1
  (center (p2 ph1) (p22 ph1):f3 )
  DELTA1
  p17:gp1
  d16

  ; t1 evolution
  if "l0 %2 == 1"
     {
     (p21 ph11):f3  ; 1 2 3 0
     }
  else
     {
     (p21 ph12):f3  ; 1 0 3 2
     }
#ifndef ONE_D
  d0 gron0
  d0 gron0*-1
  10u groff
#endif

  ; back transfer
  if "l0 %2 == 1"
     {
     (center (p1 ph1 4u p1 ph1):f1 (p21 ph1):f3 ) ; x/x
     }
  else
     {
     (center (p1 ph1 4u p1 ph2):f1 (p21 ph1):f3 ) ; x/-x
     }
  p16:gp2
  d16
  DELTA2
  if "l0 %2 == 1"
     {
     (p11:sp3 ph2):f1 ; FD
     }
  else
     {
     (p11:sp5 ph2):f1 ; stronger FD
     }
  4u pl1:f1

  (center (p2 ph1) (p22 ph1):f3 )

  4u
  (p11:sp4 ph2):f1  ; FB
  DELTA3
  p16:gp2
  d16 pl16:f3
  4u BLKGRAD

;  go=2 ph31 cpd3:f3
;  d11 do:f3 mc #0 to 2 
;     F1EA(iu0, id0)

  gosc ph31 cpd3:f3
  30u do:f3
  ; collect half of phase cycle
  ;lo to 2 times 2

  ; store data and move to next buffer
  30u st
  lo to 3 times nbl

  30u
  lo to 4 times ns

  ; save buffer contents to disk
  d11 wr #0 if #0
  30u zd
  
  ; inner loop (echo/anti-echo)
  30u iu0
  lo to 5 times 2

  ; outer loop (d0)
  30u id0
  lo to 6 times l2

exit 
  

ph1=0
ph2=2
ph3=0
ph11=1 3 2 0
ph12=1 3 0 2
;ph31=0 2 3 1  ; ZQ
ph31=0 2 1 3  ; DQ


;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;sp1: f1 channel - 90 degree water flip-down (sinc1000)
;sp2: f1 channel - 90 degree water flip-down (sinc1000)
;sp3: f1 channel - 90 degree water flip-down (sinc1000)
;sp4: f1 channel - 90 degree water flip-back (sinc1000)
;sp5: f1 channel - 90 degree water flip-down (sinc1000)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse                       [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                         [6 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d16: delay for homospoil/gradient recovery
;d26 : 1/(4J)YH
;cnst4: = J(YH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;ns: 1 * n
;nbl: 4
;ds: 16
;td1: 4
;FnMODE: echo-antiecho


;for z-only gradients:
;gpz0: 3%
;gpz1: 30%
;gpz2: 50%
;gpz3: 11%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100


;use AU-program xfshear to process data



;$Id: trosyzqgpphwg,v 1.8 2012/01/31 17:49:31 ber Exp $
