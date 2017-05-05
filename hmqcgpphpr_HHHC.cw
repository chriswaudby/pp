; For measurement of HC/HH D/D CCR
; Sun & Tugarinov (2012) [Fig. 2B]
;
; IP and AP spectra interleaved in F1
; set td1 = 2 * number of points in vdlist
; set vdlist = 0, 4, 8... ms
;
; F2 = 13C dimension
; set td2 = 1, and -DONE_D to acquire pseudo-2D

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

"p2=p1*2"
"p4=p3*2"
"d4=p4"

"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf2/2"

# ifdef ONE_D
"d0=0.1u"
# else
"d0=in0/2-0.63662*p3-p1"
# endif /*ONE_D*/

"DELTA1=2m-p16-d16-larger(p1,p3)"
"DELTA2=1m-p19-d16-larger(p1,p3)"
"DELTA3=500u-p17-d16-larger(p1,p3)"
"DELTA4=DELTA1-p4-de"

"l2=0"
"l3=td2/2"
"l6=td1/2"

aqseq 312

1 ze 
  d11
2 d11 
3 4u
4 4u
5 4u
6 4u

  "TAU=vd/2-larger(p3,p2+0.1u)"

  4u BLKGRAD

  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  d1 cw:f1 ph1
  d13 do:f1
  d12 pl1:f1 pl2:f2
  30u fq=0:f1
  50u UNBLKGRAD
  (p3 ph1):f2    ; crush eq'm 13C magnetisation
  d13
  p16:gp1
  d16*2 

  (p1 ph1):f1
  p16:gp2
  d16
  DELTA1
  (center (p2 ph1):f1 (p4 ph1):f2 )
  DELTA1
  p16:gp2
  d16

  (p3 ph11):f2
  p19:gp3
  d16
  DELTA2
  (center (p2 ph1):f1 (p4 ph1):f2 )
  DELTA2
  p19:gp3
  d16

  (p3 ph12):f2
  d0
  (p2 ph1):f1
  d0

  (lalign (p1 ph1):f1 (p3 ph13):f2 ) 
  p19:gp4
  d16
  DELTA2
  (center (p2 ph1):f1 (p4 ph1):f2 )
  DELTA2
  p19:gp4
  d16

  (p3 ph14):f2
  0.1u
  (p3 ph15):f2
  p17:gp5
  d16
  DELTA3
  (center (p2 ph1):f1 (p4 ph1):f2 )
  DELTA3
  p17:gp5
  d16

  (p3 ph16):f2
  TAU       ;relaxation delay
  (center (p1 ph1 0.1u p2 ph2 0.1u p1 ph1):f1 (p4 ph1):f2 )
  TAU

  (p1 ph1)

if "l2 % 2 == 0"
{
  p16:gp6
  d16
  DELTA1
  (center (p2 ph1):f1 (p4 ph1):f2 )
  DELTA4
  p16:gp6
  d16
  d4 BLKGRAD
}
else
{
  p16:gp6
  d16
  DELTA1
  (center (p2 ph1):f1 (d4) )
  DELTA4
  p16:gp6
  d16 BLKGRAD
  (p4 ph1):f2
}

  go=2 ph31 

  d11 wr #0 if #0 zd

  2u iu2
  2u ip31
  lo to 3 times 2

  2u ivd
  2u rp31
  lo to 4 times l6

# ifndef ONE_D
  4u ru2
  4u ip13
  4u ip14
  lo to 5 times 2

  4u id0
  lo to 6 times l3
# endif /*ONE_D*/

;  d11 mc #0 to 2
;      F1QF(ivd)
;      F2PH(ip13 & ip14, id0)

exit 
  
  
ph1= 0 
ph2= 1 
ph11=0 2
ph12=1 1 3 3
ph13=1 1 1 1 3 3 3 3
ph14=1 1 1 1 1 1 1 1
     3 3 3 3 3 3 3 3
ph15=1
ph16=(8) 1 1 5 5      ; for L1 selection
;ph17=(8) 3 3 7 7     ; for L3 selection
ph31=0 2 2 0 0 2 2 0
     2 0 0 2 2 0 0 2


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse [1000 usec]
;p17: gradient pulse [300 usec]
;p19: gradient pulse [50 usec]
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery          [200 usec]
;cnst21: frequency in Hz for off-res presat
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 16 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ

;for z-only gradients:
;gpz1: 31%
;gpz2: 13%
;gpz3: 11%
;gpz4: 7%
;gpz5: 29%
;gpz6: 17%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SINE.10
;gpnam4: SINE.10
;gpnam5: SINE.10
;gpnam6: SINE.10

