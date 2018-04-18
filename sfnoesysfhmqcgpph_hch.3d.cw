;3D H[N/M]CMHM SF-NOESY-SFHMQC
;for amide/methyl-methyl NOES
; Rossi 2016 JBNMR

;F1(H,t1) ---NOE--> F1(H) -> F2(C[mq],t2) -> F1(H,t3)
;
;Indirect evolution order is t2, t1 (13C, 1Hnoe)
;Delays adjusted for zero first-order phase correction in acqusition dim


;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d11=45m"
"d12=20u"
"d13=4u"

;------------options for first (in transfer pathway) 1H (F1)
"d0=3u"
"in0=inf1"

;------------options for second (in transfer pathway) 13C (F2)
"in10=inf2"
"d10=in10/2-p3*4/3.1415"


"spoff13=bf1*(cnst19/1000000)-o1" ; H[N] fd (p29)
"spoff14=bf1*(cnst19/1000000)-o1" ; H[N] 180 (p30)
"spoff15=bf1*(cnst19/1000000)-o1" ; H[N] fb (p29)
"spoal13=1"
"spoal14=0.5"
"spoal15=0"

"spoff23=bf1*(cnst20/1000000)-o1" ; H[M] fd (p39)
"spoff24=bf1*(cnst20/1000000)-o1" ; H[M] 180 (p40)
"spoal23=1"
"spoal24=0.5"

"TAU=d8-p16*2-d16*2-p3-14u"


"DELTA1=d2-p16-d16"
"DELTA2=p39*cnst39-4u-de"
"acqt0=de"

aqseq 321	; for info only


1 ze
  d11 pl12:f2 pl16:f3
2 d11 do:f2 do:f3
  4u BLKGRAD

  d1 
  50u UNBLKGRAD
  p16:gp1
  d16

;-------------------------start 1H evolution

  20u pl16:f3
  4u cpd3:f3 ; 15N cpd
  (p29:sp13 ph11):f1
  d0
  (p30:sp14 ph13):f1
  3u
  (p29:sp15 ph14):f1
  4u do:f3

  ; purge residual magnetisation
  p16:gp1
  d16 pl2:f2
  (p3 ph1):f2
  10u
  p16:gp1*0.71
  d16


  TAU


;------------------------start second 13C HMQC element
  (p39:sp23 ph1):f1
  p16:gp2
  d16
  (center (p40:sp24 ph2):f1 (DELTA1 p3 ph12 d10 p3 ph15 DELTA1):f2 )
  DELTA2
  p16:gp2
  d16 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2

  d11 do:f2 mc #0 to 2
	F2PH(ip15, id10)
	F1PH(rp15 & rd10 & ip14, id0)

  4u BLKGRAD

exit


ph1= 0
ph2= 1
ph11=0 2
ph12=0 0 2 2
ph13=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph14=0
ph15=0
ph31=0 2 2 0 2 0 0 2


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;pl16: f3 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p29: f1 channel - 90 degree shaped pulse for excitation (amide)
;                      Pc9_4_90.1000 (90o)    (2657 us at 700 MHz)
;p30: f1 channel - 180 degree shaped pulse for refocussing (amide)
;                      Reburp.1000               (1943 us at 700 MHz)
;p39: f1 channel - 120 degree shaped pulse for excitation (methyl)
;                      Pc9_4_120.1000 (120o)    (2314 us at 700 MHz)
;p40: f1 channel - 180 degree shaped pulse for refocussing (methyl)
;                      Reburp.1000               (1693 us at 700 MHz)
;sp13: f1 channel - shaped pulse 90 degree (amide)
;                   (Pc9_4_90.1000 )
;sp14: f1 channel - shaped pulse 180 degree (Reburp.1000) (amide)
;sp23: f1 channel - shaped pulse 120 degree (methyl)
;                   (Pc9_4_120.1000 )
;sp24: f1 channel - shaped pulse 180 degree (Reburp.1000) (methyl)
;sp25: f1 channel - shaped pulse 120 degree (methyl)
;                   (Pc9_4_120.1000 )
;cnst19: H(N) chemical shift (offset, in ppm) [8.2 ppm]
;cnst20: H(methyl) chemical shift (offset, in ppm) [0.7 ppm]
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529

;p16: homospoil/gradient pulse                       [1 msec]
;p22: f3 channel - 180 degree high power pulse
;p28: f1 channel - trim pulse                        [1 msec]
;d0 : incremented delay (4D)
;d10: incremented delay (4D)
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d8 : mixing time
;d11: delay for disk I/O                             [45 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(CH)
;inf1: 1/SW(H)
;inf2: 1/SW(C) = DW(C)
;in0: 1/(SW(H)) = DW(C)
;in10: 1/(2 * SW(C)) = 0.5 * DW(C)
;NS: 8 * n
;DS: 16
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI (or States) in F1
;FnMODE: States-TPPI (or States) in F2
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence
;cpd3: decoupling according to sequence defined by cpdprg2
;pcpd3: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz2: 29%
;gpz3: 23%
;gpz4: 13%
;gpz5: 43%
;gpz6: 19%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SINE.50
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam5: SINE.50
;gpnam6: SINE.100

                                          ;preprocessor-flags-start
;LABEL_CN: for 15N decoupling during indirect 13C evolution periods
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
;NOE_SAT: for water saturation during NOE mixing time
;F2_plane: for zero 13C phase evolution in F3
;F3_plane: for zero 13C phase evolution in F2
;F2_SINGLEDWELL: for single-dwell first-point delay in F2
;F3_SINGLEDWELL: for single-dwell first-point delay in F3
;TRIMPULSE: to apply trim-pulses in second 13C HMQC element, set p28
;NUS: for non-uniform sampling (Topspin 3)
                                          ;preprocessor-flags-end
