; multiplet-resolved 13C T1p measurement
;   based on hsqcphpr_1Hcoupled.3.cw
; run as pseudo-3D
; set td1 = 4 * number of relaxation times
;
; to separate into sub-components:
;   aaa = A + 2B + 2C + D
;   aab = A +  B -  C - D
;   abb = A -  B -  C + D
;   bbb = A - 2B + 2C - D 
;
; May 2017
;
; with off-resonance presat
; ZZ/crusher periods, clean-up gradient pairs
; (90,-180) phase correction
; use baseopt
;
;hsqcphpr
;avance-version (07/04/04)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive
;with decoupling during acquisition
;
;G. Bodenhausen & D.J. Ruben, Chem. Phys. Lett. 69, 185 (1980)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Delay.incl>
#include <Grad.incl>

define list<delay> Trelax = <$VDLIST>

"p2=p1*2"
"d2=p2"
"p4=p3*2"
"p22=p21*2"
"d4=1s/(cnst2*4)"
"d11=30m"
"d12=20u"
"d13=4u"


"TAU1=0.1u"
"TAU2=1s/(cnst2*12)-0.6366*p3-larger(p1,p3)"
"TAU3=1s/(cnst2*6)-0.6366*p3-larger(p1,p3)"
"TAU4=1s/(cnst2*4)-0.6366*p3-larger(p1,p3)"

"in0=inf2"
"d0=in0/2-2*0.6366*p3"

"DELTA=d4-p16-d16-larger(p1,p3)-0.6366*p1"
"DELTA1=d4-p19-d16-p10-p1-4u-0.6366*p1"
"DELTA2=d4-p19-d16-p10-p1-12u"
"acqt0=0"

; calculate offset for WFB
"spoff1=cnst21-o1"

; loop counter for interleaving
"l0=0"

; loop counters for acquisition
"l1 = td1 / 4"
"l2 = td2 / 2"

aqseq 312

1 ze
  d11 pl12:f2
2 d11 do:f2
3 4u
4 4u
5 4u
6 4u

  ; off-resonance presat
  30u pl9:f1
  30u fq=cnst21(bf hz):f1
  d1 cw:f1 ph1
  30u do:f1
  30u fq=0:f1

  ; purge equilibrium 13C
  30u UNBLKGRAD
  4u pl1:f1 pl2:f2
  (p3 ph1):f2
  p16:gp0
  d16

  ; begin main sequence
  (p1 ph1)
  p16:gp1
  d16
  DELTA
  (center (p2 ph1) (p4 ph1):f2 )
  DELTA
  p16:gp1
  d16
  (p1 ph11)

  ; zz purge
  p16:gp2
  d16

  ; multiplet filter
  if "l0 % 4 == 0"
  {
  (p3 ph1):f2
  (center (p1 ph2 p2 ph1 p1 ph2):f1 (TAU1 p4 ph1 TAU1):f2 )
  (p3 ph1):f2
  }

  if "l0 % 4 == 1"
  {
  (p3 ph1):f2
  (center (p1 ph2 p2 ph1 p1 ph2):f1 (TAU2 p4 ph1 TAU2):f2 )
  (p3 ph2):f2
  }

  if "l0 % 4 == 2"
  {
  (p3 ph1):f2
  (center (p1 ph2 p2 ph1 p1 ph2):f1 (TAU3 p4 ph1 TAU3):f2 )
  (p3 ph3):f2
  }

  if "l0 % 4 == 3"
  {
  (p3 ph1):f2
  (center (p1 ph2 p2 ph1 p1 ph2):f1 (TAU4 p4 ph1 TAU4):f2 )
  (p3 ph4):f2
  }

  ; 13C spin-lock and t1
  (p3 ph13):f2
  0.1u pl8:f2
  Trelax cw:f2 ph2
  0.1u do:f2
  d0 pl2:f2
  (p3 ph14):f2

  ; zz purge
  p16:gp3
  d16

  ; final inept
  (p1 ph1)
  p19:gp4
  d16
  DELTA1
  (p10:sp1 ph3):f1
  4u pl1:f1
  (center (p2 ph1) (p4 ph1):f2 )
  4u
  (p10:sp1 ph3):f1
  DELTA2
  p19:gp4
  d16
  4u BLKGRAD
  4u pl12:f2

  go=2 ph31 cpd2:f2
  d11 do:f2 wr #0 if #0 zd

  4u Trelax.inc
  lo to 3 times l1

  4u iu0
  lo to 4 times 4

  4u ip13
  lo to 5 times 2

  4u id0
  lo to 6 times l2


exit


ph1=0
ph2=1
ph3=2
ph4=3
ph11=1 1 1 1 3 3 3 3
ph13=0 2
ph14=0 0 2 2
ph31=0 2 2 0 2 0 0 2


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p10 : f1 channel - 90 degree selective pulse [1000 usec]
;sp1 : f1 channel - 90 degree WFB (p10)
;d0 : incremented delay (2D)
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;cnst2: = J(XH)
;cnst21: off-resonance presaturation frequency (bf hz)
;inf1: 1/SW(X) = DW(X)
;in0: 1/SW(X) = DW(X)
;nd0: 1
;NS: 8 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence

;for z-only gradients:
;gpz0: 46 %
;gpz1: 13 %
;gpz2: 17 %
;gpz3: 33 %
;gpz4: 29 %

;gradients
;p16: 1000u
;p19: 300u

;use gradient files:
;gpnam0: SMSQ10.100
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SINE.10
