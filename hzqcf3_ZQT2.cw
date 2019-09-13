; Gradient-selected HZQC (gsHZQC) methyl-TROSY with  methyl filter and ZQ Hahn-echo relaxation delay
; For recording 1H-13C ZQ relaxation rates in fully protonated methyl groups
; With or without 13C polarization enhancement
; Reference: Gill, ML and Palmer, III, AG, Journal of Biomolecular NMR (2011) XX:XXX-XXX, Figure 3
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=

prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

; Ensure the -DNOPOL is set if ZGOPTS is empty
#ifndef POL1
#ifndef POL2
#ifndef NOPOL
#define NOPOL
#endif
#endif
#endif

; Unset other flags if using one of polarization schemes
#ifdef POL2
#undef POL1
#undef NOPOL
#endif

#ifdef POL1
#undef POL2
#undef NOPOL
#endif


"cnst21=o1/bf1" ; acquisition frequency should be set the same as O1P
"p11=p1*2"      ; 1H 180degree pulse at hard power
"p21=p2*2"      ; 13C 180degree pulse at hard power
"d11=10u"       ; hardware delay
"d16=250u"      ; gradient recovery delay
"d28=p21"       ; compensation for 13C 180
"d29=p1*4+8u"   ; compensation for 1H composite pulse
"TAU=3.91m"     ; 1/(2JCH), JCH=128 Hz
"TAU2=0.75m"    ; sin(2piJCH*TAU2) = 3^(–1/2)
"in0=inf1"      ; t1 increment
"d0=in0/2"      ; initial t1 delay (1/2 dwell)
"l0=td1/4"      ; number of t1 complex points
"l1=0"          ; initialize loop counter

"DELTA1=TAU/4"
"DELTA2=(TAU/4)-p16-d16-d11"
"DELTA7=TAU*cnst2/2"
"DELTA8=(TAU*cnst2/2)-p16-d16-d11"

#ifdef NOPOL
"DELTA3=TAU-d11-p1*0.6366"
"DELTA4=TAU-p16-d16-d11"
#endif

#ifdef POL1
"DELTA3=TAU/2-1*TAU2-larger(p1,p2)-p1*0.6366-p17-d16"
"DELTA4=TAU/2+1*TAU2-larger(p1,p2)"
"DELTA5=TAU/2+1*TAU2-larger(p1,p2)-p17-d16-4u"
"DELTA6=TAU/2-1*TAU2-d16-p16-d11-4u-larger(p1,p2)"
#endif

#ifdef POL2
"DELTA3=TAU-p1*0.6366"
"DELTA4=TAU-p16-d16-d11"
"DELTA5=TAU/2+TAU2-p1-p17-d16"
"DELTA6=TAU/2-TAU2-p1-p17-d16-4u"
#endif


1 ze
  d11 BLKGRAD
2 30m
3 d11 do:f2 fq=cnst20(bf ppm):f1 
  d11 pl9:f1
  d1 cw:f1 ph0 pl2:f2
  4u do:f1
  d11 pl1:f1 fq=cnst21(bf ppm):f1
  d11 UNBLKGRAD
  
#ifdef NOPOL
  ;***purge***
  (p2 ph0):f2
  4u
  p15:gp1
  d16
  ;***evolve AP to MQ***
  (p1 ph1):f1
  d11 BLKGRAD
  DELTA3
  if "l1 % 4 == 0" goto 4
  if "l1 % 4 == 1" goto 5
  if "l1 % 4 == 2" goto 6
  if "l1 % 4 == 3" goto 7
#endif

#ifdef POL1
  p15:gp1
  d16
  ;***transfer polarization from 13C to 1H***
  (p2 ph2):f2
  TAU2
  TAU2
  (p1 ph1):f1
  ;***evolve AP to MQ***
  DELTA3
  p17:gp6
  d16
  (center (p11 ph20):f1 (p21 ph2):f2 )
  4u
  p17:gp6
  d16
  DELTA5
  if "l1 % 4 == 0" goto 6
  if "l1 % 4 == 1" goto 7
  if "l1 % 4 == 2" goto 4
  if "l1 % 4 == 3" goto 5
#endif

#ifdef POL2
  p15:gp1
  d16
  ;***transfer polarization from 13C to 1H***
  (p2 ph2):f2
  DELTA5
  p17:gp6
  d16
  (p11 ph20):f1
  4u
  p17:gp7
  d16
  DELTA6
  (ralign (p1 ph1):f1 (p21 ph7):f2 )
  ;***evolve AP to MQ***
  DELTA3
  if "l1 % 4 == 0" goto 4
  if "l1 % 4 == 1" goto 5
  if "l1 % 4 == 2" goto 6
  if "l1 % 4 == 3" goto 7
#endif

 4  d29
    ;***create MQ magnetization***
    (p2 ph2):f2
    4u
    (p21 ph20):f2
    ;***J-filter***
    DELTA1
    4u
    ;***Relaxation delay (T/2)***
    DELTA7
    ;***t1***
    d0
    (p1 ph10 4u p11 ph11 4u p1 ph10):f1
    (p21 ph20):f2
    ;***J-filter***
    DELTA1
    d29
    ;***Relaxation delay (T/2)***
    DELTA8
    d11 UNBLKGRAD
    p16:gp2 ; ZQ encode gradient
    d16
    ;***create AP magnetization***
    (p2 ph5):f2
    8u
 goto 8

 5  d29
    ;***create MQ magnetization***
    (p2 ph2):f2
    4u
    ;***J-filter***
    DELTA1
    (p21 ph20):f2
    4u
    ;***Relaxation delay (T/2)***
    DELTA7
    ;***t1***
    d0
    (p1 ph10 4u p11 ph11 4u p1 ph10):f1
    ;***J-filter***
    DELTA1
    (p21 ph20):f2
    d29
    ;***Relaxation delay (T/2)***
    DELTA8
    d11 UNBLKGRAD
    p16:gp2 ; ZQ encode gradient
    d16
    ;***create AP magnetization***
    (p2 ph5):f2
    8u
 goto 8

 6  d29
    ;***create MQ magnetization***
    (p2 ph2):f2
    4u
    (p21 ph20):f2
    ;***Relaxation delay (T/2)***
    DELTA7
    4u
    ;***J-filter***
    DELTA1
    (p1 ph10 4u p11 ph11 4u p1 ph10):f1
    (p21 ph20):f2
    d0
    ;***Relaxation delay (T/2)***
    DELTA7
    d29
    ;***J-filter***
    DELTA2
    d11 UNBLKGRAD
    p16:gp2 ; ZQ encode gradient
    d16
    ;***create AP magnetization***
    (p2 ph5):f2
    8u
 goto 8

 7  d29
    ;***create MQ magnetization***
    (p2 ph7):f2
    4u
    ;***Relaxation delay (T/2)***
    DELTA7
    (p21 ph20):f2
    4u
    ;***J-filter***
    DELTA1
    (p1 ph10 4u p11 ph11 4u p1 ph10):f1
    ;***t1***
    d0
    ;***Relaxation delay (T/2)***
    DELTA7
    (p21 ph20):f2
    d29
    ;***J-filter***
    DELTA2
    d11 UNBLKGRAD
    p16:gp3 ; DQ encode gradient
    d16
    ;***create AP magnetization***
    (p2 ph6):f2
    8u
 goto 8

#ifdef POL1
 8  DELTA4
    (center (p11 ph10):f1 (p21 ph20):f2 )
    4u
    p16:gp5 ; SQ decode gradient
    d16 pl12:f2
    d11 BLKGRAD
    DELTA6
#else
 8  p16:gp4 ; SQ decode gradient
    d16 pl12:f2
    d11 BLKGRAD
    DELTA4 
#endif

    ;***detect***
    go=2 ph30 cpd2:f2
    10u do:f2

    if "l1 % 2 == 1" goto 10

    30m
    goto 11

10  30m wr #0 if #0 zd
11  30m iu1
    lo to 2 times 4
    30m id0
    lo to 3 times l0
exit


ph0=0
ph1=0 0
ph2=0 2
ph5=0 0
ph6=2 2
ph7=2 0
ph10=0 0
ph11=1 1
ph20=0 0

#ifdef POL2
ph30=2 0
#else
ph30=0 2
#endif

;*** VARIABLES SET BY USER ***
;d1 = recycle delay
;td1 = 4 * number of complex points (e.g. for 250 t1 points, td1=1000)
;o1p = set offset to just outside of methyl region [2.0 ppm]

;p1 = 1H hard 90 @ pl1
;p2 = 13C hard 90 at pl2
;p11 = 1H hard 180 @ pl1
;p21 = 13C hard 180 @ pl2
;p15 = pulse for gz1 [1ms]
;p16 = pulse for gz2-gz5 [500us]
;p17 = pulse for gz6-gz7 [500us]

;pl0 = 120 dB
;pl1 = 1H hard power
;pl2 = 13C hard power
;pl12 = 13C broadband decoupling power
;p15 = pulse for gpz1 [1ms]
;p16 = pulse for gpz2-gpz5 [500us]
;p17 = pulse for gpz6-gpz7 [500us]

;cnst2 = number of 1/(2JCH) intervals for relaxation delay (T)
;cnst20 = water offset [4.7ppm]

;gpz1 7.5 G/cm (15%)
;gpz2 30 G/cm (60%)
;gpz3 18 G/cm (36%)
;gpz4 -22.5 G/cm (-45%)
;gpz5 22.5 G/cm (45%)
;gpz6 5 G/cm (10%)
;gpz7 -5 G/cm (-10%)

;gpnam1: SINE.100
;gpnam2: SINE.50
;gpnam3: SINE.50
;gpnam4: SINE.50
;gpnam5: SINE.50
;gpnam6: SINE.50
;gpnam7: SINE.50

; F1 frequency discrimination has been coded for echo/antiecho, thus
; FnMODE: undefined

; preprocessor flags, use ONLY ONE of the following:
; NOPOL : to purge 13C polarization, this is the default if nothing is set
;         option -DNOPOL (eda: ZGOPTNS)
; POL1  : to use 13C polarization enhancement scheme in Figure 1B
;         option -DPOL1 (eda: ZGOPTNS)
; POL2  : to use 13C polarization enhancement scheme in Figure 1C
;         option -DPOL2 (eda: ZGOPTNS)

; Data processing: For each consecutive block of four data rows,
;                  add the first to the second and the third to the fourth, 
;                  then process using Rance-Kay (Echo-Antiecho) for frequency discrimination
