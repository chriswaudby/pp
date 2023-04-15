/* 13CH3 STE measurement using DQ'/QQ coherences (overall TQ encoding)
  based on hsqc_1D_3Q_diff_600_lek_cp
  Used to record 1D diffusion using 13CH3 methyl by creating 1H 3Q
    Magnetization is of the form and IZCZ during diff time
  Written by LEK on Jan 29, 2016
  72 step phase cycle
  Modified by LEK on April 9, 2017 to include the option of having magnetization   of
the form Iz and IzIzIz during diffusion delay
  Modified by LEK on April 20, 2017 to include water gate (-Dw_gate)
*/

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

;Define phases
#define zero ph=0.0
#define one ph=90.0
#define two ph=180.0
#define three ph=270.0

;Define Pulses
define pulse pwh
"pwh=p1"  ; 1H hard pulse at power pl1

define pulse pwc
"pwc=p3"  ; 13C hard pulse at power pl2

;Define delays
"in0=inf1/2"
"d11=30m"

define list<gradient> diff=<Difframp>

define delay taua
"taua=d5"   ; < 1/(4JCH) 1.8 ms

define delay taub
"taub=d6"   ; = 1/(4JCH) 2.0 ms

define delay T_diff
"T_diff = d20"

define delay hscuba
"hscuba=30m"

/* Assign cnsts to check validity of parameter changes */
#ifdef fsat
  "cnst10 = plw10" ; tsatpwr - set to max of 0.00005W
#endif

"cnst21=plw21" ; dpwr pl21 - set max at 2.5W

#ifndef w_gate
  "acqt0=-pwh*2.0/PI - 2u"  ; select 'DIGIMOD = baseopt' to execute
#else
  "acqt0=-2u"  ; select 'DIGIMOD = baseopt' to execute
#endif

/*************************************/
/* Define macro for 3-9-19 WATERGATE */
/*************************************/
#ifdef w_gate
#define Watergate (center (pwh*0.231 ph26 d19*2 pwh*0.692 ph28 d19*2 pwh*1.462 ph26 d1 9*2 pwh*1.462 ph26 d19*2 pwh*0.692 ph28 d19*2 pwh*0.231 ph26):f1 (pwc*2 ph26):f2)
#endif

1 ze

; check validity of parameters

#ifdef fsat
  if "cnst10 > 0.001" {
    2u
    print "error: presat power pl10 is too large"
    goto HaltAcqu
  }
#endif
if "ns % 72 !=0" {
  2u
  print "error: ns must be a multiple of 72"
  goto HaltAcqu
}
if " cnst21 > 2.5" {
  2u
  print "error: dpwr pl21 too large"
  goto HaltAcqu
}
if "pwc > 20u" {
  2u
  print "error: pwc too large < 20 us"
  goto HaltAcqu
}

2 d11 do:f2

#ifdef scramble
  ; scramble any residual 1H magnetization
  20u UNBLKGRAD   ; dly 20u, unblank gradients and lock hold
  2u pl1:f1
  (pwh one):f1
  2u
  p56:gp6*0.25    ; gradient 6 * 0.25
  d16
  2u pl1:f1
  (pwh zero):f1
  2u
  p56:gp6*0.5     ; gradient 6*0.5
  d16
  20u BLKGRAD     ; dly 20u, blank gradients
#endif

#ifdef fsat
  ; zgoptn -Dfsat
  20u fq=cnst1:f1         ; jump from o1 to water
  4u pl10:f1              ; power(tsatpwr) for presaturation
  d1 cw:f1 zero           ; Hcw(d1)x
  4u do:f1                ; cw off
  2u pl1:f1               ; power(tpwr)
#ifdef fscuba
    /* Scuba pulse sequence */ 
    hscuba                ; delay(hscuba)
    (pwh zero):f1         ; H 90x180y90x
    (pwh*2 one):f1
    (pwh zero):f1
    hscuba                ; delay(hscuba)
#endif  /* fscuba */
#else   /* if fsat is no */
  2u pl1:f1               ; power(tpwr)
  d1                      ; delay(d1)
#endif    /* end if fsat */

2u pl1:f1
2u pl2:f2
2u pl31:f3

20u UNBLKGRAD   ; dly 20u, unblank gradients and lock on hold
10u fq=0:f1     ; jump back to o1

(pwc zero):f2   ; C90x - To destroy 13C Boltzman magnetization
  2u
  (p50:gp0)     ; gradient 0
  d16


; This is the real start

(pwh ph1):f1                    ; H90ph1

2u
p51:gp1                         ; gradient 1
d16

; subtract effects of both 1H 90
"DELTA = taua - 2.0u - p51 - d16 - pwc - pwh*4/PI - 0.2u"
DELTA                                    ; delay 1/4JCH

(center (pwh*2 ph1):f1 (pwc*2 ph26):f2) ; H180x,C180x

"DELTA = taua - 2.0u - p51 - d16 - pwc - 2u"
DELTA                                    ; delay 1/4JCH

2u
p51:gp1           ; gradient 1
d16

2u pl2:f2
(pwc ph3):f2      ; C90ph3

2u
p61:gp11          ; gradient 11
d16

"DELTA = taub - 2.0u - p61 - d16"
DELTA             ; delay 1/4JCH

(center(pwh*2 ph1):f1 (pwc ph27 pwc*2 ph26 pwc ph27):f2) ; H180x,C180x

"DELTA = taub - 2.0u - p61 - d16"
DELTA             ; delay 1/4JCH

2u
p61:gp11          ; gradient 11
d16

(pwh ph1):f1      ; H90t2 - FIRST DIFFUSION ENCODING PERIOD

2u
p53:gp3*diff      ; gradient 3 diffusion encoding gradient
d17

(pwc ph27 pwc*2 ph26 pwc ph27):f2   ; C180x

2u
p53:gp3*diff      ; gradient 3 diffusion encoding gradient
d17

(pwh ph27 pwh*2 ph26 pwh ph27):f1

2u
p53:gp3*-1.0*diff ; -gradient 3 diffusion encoding gradient
d17

(pwc ph27 pwc*2 ph26 pwc ph27):f2   ; C180x

2u
p53:gp3*-1.0*diff ; -gradient 3 diffusion encoding gradient
d17


#ifndef Iz
  ; magnetization diffuses as 2IzCz
  (pwh zero):f1          ; H90x

  2u
  p54:gp4               ; gradient 4
  d16

  "DELTA = taub - 2.0u - p54 - d16"
  DELTA                 ; delay 1/4JCH

  (center(pwh ph27 pwh*2 ph26 pwh ph27):f1 (pwc ph27 pwc*2 ph26 pwc ph27):f2) ; H180x,C180x

  "DELTA = taub - 2.0u - p54 - d16"
  DELTA                 ; delay 1/4JCH

  2u
  p54:gp4               ; gradient 4
  d16

  (pwc zero):f2         ; C90x
  1.0u
  (pwh one):f1          ; H90y
#else
  ; magnetisation diffuses as Iz
  (pwh ph3):f1          ; H90ph3
  
  2u
  p54:gp4               ; gradient 4
  d16

  "DELTA = taub - 2.0u - p54 - d16"
  DELTA                 ; delay 1/4JCH

  (center(pwh ph27 pwh*2 ph26 pwh ph27):f1 (pwc ph27 pwc*2 ph26 pwc ph27):f2) ; H180x,C180x

  "DELTA = taub - 2.0u - p54 - d16"
  DELTA                 ; delay 1/4JCH
  
  2u
  p54:gp4               ; gradient 4
  d16

  (pwc zero):f2         ; C90x

  2u
  p54:gp4*0.5           ; gradient 4 * 0.5
  d16

  "DELTA = taub - 2.0u - p54 - d16"
  DELTA                 ; delay 1/4JCH

  (center(pwh ph27 pwh*2 ph26 pwh ph27):f1 (pwc ph27 pwc*2 ph26 pwc ph27):f2) ; H180x,C180x
  
  "DELTA = taub - 2.0u - p54 - d16 - pwh*2.0/PI"
  DELTA                 ; delay 1/4JCH

  2u
  p54:gp4*0.5           ; gradient 4 * 0.5
  d16

  (pwh zero):f1
#endif

#ifndef Iz
  ; magnetization diffuses as 2IzCz
  "DELTA = T_diff - taub*4 - 2u - p56 - d16"
#else
  "DELTA = T_diff - taub*8 - 2u - p56 - d16"
#endif

#ifndef fsat_diff
  DELTA
#endif

#ifdef fsat_diff
  ; zgoptn -Dfsat
  10u fq=cnst1:f1       ; jump from o1 to water
  2u pl10:f1            ; power(tsatpwr) for presaturation
  DELTA cw:f1 zero      ; Hcw(d1)x
  2u do:f1              ; cw off
  2u pl1:f1             ; power(tpwr)
  10u fq=0:f1           ; jump back to i1
#endif

2u
p56:gp6                 ; gradient 6
d16

#ifndef Iz
  (pwh one):f1          ; H90y
  1.0u
  (pwc zero):f2          ; C90x

  2u
  p55:gp5               ; gradient 5
  d16

  "DELTA = taub - 2.0u - p55 - d16"
  DELTA                 ; delay 1/4JCH

  (center(pwh*2 zero):f1 (pwc ph27 pwc*2 ph26 pwc ph27):f2) ; H180x,C180x

  "DELTA = taub - 2.0u - p55 - d16"
  DELTA                 ; delay 1/4JCH

  2u
  p55:gp5               ; gradient 5
  d16

  (pwh zero):f1          ; H90ph6
#else
  (pwh ph14):f1         ; H90ph14

  2u
  p55:gp5*0.5           ; gradient 5 * 0.5
  d16
  
  "DELTA = taub - 2.0u - p55 - d16 - pwh*2.0/PI - 1u"
  DELTA                 ; delay 1/4JCH

  (center(pwh*2 ph14):f1 (pwc ph27 pwc*2 ph26 pwc ph27):f2) ; H180x,C180x
  
  "DELTA = taub - 2.0u - p55 - d16"
  DELTA                 ; delay 1/4JCH

  2u
  p55:gp5*0.5           ; gradient 5 * 0.5
  d16

  (pwc zero):f2         ; C90x

  2u
  p55:gp5               ; gradient 5
  d16

  "DELTA = taub - 2.0u - p55 - d16"
  DELTA

  (center(pwh*2 ph14):f1 (pwc ph27 pwc*2 ph26 pwc ph27):f2) ; H180x,C180x

  "DELTA = taub - 2.0u - p55 - d16"
  DELTA
  
  2u
  p55:gp5               ; gradient 5
  d16

  (pwh ph14):f1         ; H90ph14
#endif

2u
p53:gp3*diff      ; gradient 3 diffusion encoding gradient
d17

(pwc ph27 pwc*2 ph26 pwc ph27):f2   ; C180x

2u
p53:gp3*diff      ; gradient 3 diffusion encoding gradient
d17

(pwh ph27 pwh*2 ph26 pwh ph27):f1    ; H180 composite

2u
p53:gp3*diff      ; gradient 3 diffusion encoding gradient
d17

(pwc ph27 pwc*2 ph26 pwc ph27):f2   ; C180x

2u
p53:gp3*diff      ; gradient 3 diffusion encoding gradient
d17

(pwh ph2):f1          ; H90ph2

2u
p58:gp8               ; gradient 8
d16

"DELTA = taub - 2.0u - p58 - d16"
DELTA

(center (pwh ph27 pwh*2 ph26 pwh ph27):f1 (pwc ph27 pwc*2 ph26 pwc ph27):f2) ; H180x,C180x

"DELTA = taub - 2.0u - p58 - d16"
DELTA

2u
p58:gp8                ; gradient 8
d16

#ifndef w_gate
  (pwc zero):f2        ; C90x

  2u
  p60:gp10             ; gradient 10
  d16
  
  "DELTA = taua - 2.0u - p60 - d16 - pwc + pwh*2/PI"
  DELTA   ; delay 1/4JCH
  
  (center(pwh*2 ph26):f1 (pwc*2 ph26):f2) ; H180x,C180x
  
  2u pl2:f2
  "DELTA = taua - 2u - 2u - p60 - d16 - pwc - 4u" ; - pwc - pwc - pwh*2/PI"
  DELTA                                    ; delay 1/4JCH
  
  2u
  p60:gp10
  d16

  4u BLKGRAD
  ;(pwc zero):f2
  ;(pwc ph4):f2

  ;(pwh one):f1
#else
  ; w_gate is used
  (pwc zero):f2

  2u
  p60:gp10
  d16

  "DELTA = taua - 2.0u - p60 - d16 - d19*3.367 - 2u"
  DELTA
  2u

  Watergate

  "DELTA = taua - 2.0u - p60 - d16 - d19*3.367 - 4u - pwc - pwc"
  DELTA

  2u
  p60:gp10
  d16

  4u BLKGRAD
  (pwc zero):f2
  (pwc ph4):f2
#endif

2u pl21:f2

go=2 ph31 cpds2:f2

d11 do:f2 mc #0 to 2 F0(zd) ; write fid to disk
  F1QF(calgrad(diff))

HaltAcqu, 1m
exit


ph1= (6) 0 1 2 3 4 5  ; 0 60 120 180 240 300
ph2= (6) 0*12 1*12 2*12 3*12 4*12 5*12
ph3= 0*6 2*6
ph4=(12) 3 3 3 3 3 3 5  5  5  5  5  5  7 7 7 7 7 7
         9 9 9 9 9 9 11 11 11 11 11 11 1 1 1 1 1 1
ph5 = 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
ph6= (6) 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2
         3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5
ph7= 0
ph14= (6) 0 0 0 0 0 0 1 1 1 1 1 1 2 2 2 2 2 2
          3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5
ph31= (6) 0 3 0 3 0 3  3 0 3 0 3 0
          4 1 4 1 4 1  1 4 1 4 1 4
          2 5 2 5 2 5  5 2 5 2 5 2
ph26=0
ph27=1
ph28=2
ph29=3


;d1 : repetition delay
;d5 : taua ~1/4JCH
;d6 : taub = 1/4JCH 2 ms
;d11 : delay for disk i/o, 30ms
;d16 : gradient recovery delay, 200us
;d19 : delay for binomial water suppression = |1/(4*dis)| dis = Hz between o1 and water
;      d19 = (1/(4*|cnst1|)), cnst1 = distance between methyl and water (in Hz)
;d20: T_diff total diffusion time
;pl1 : tpwr - power level for pwh
;pl2 : dhpwr - power level for hard 13C pulse pwc
;pl21 : dpwr - power level for  13C decoupling cpd2
;pl31 : dpwr2 - power level for 15N cpd3
;p1 : pwh
;p3 : pwc
;p31 : pwn at dpwr2 for 15N decoupling during 2*TC(~29ms)
;cpd2 : 13C decoupling according to program defined by cpdprg2
;pcpd2:  13C 90 degree pulse at pl21 for cpd2
;spnam23 : file name for chirp pulse during inepts
;p53 : duration of diffusion encoding gradient pulse
;GPZ3 : set to 100 % for gradient encoding
;cnst1 : water - 01 (Hz)
;cnst12 : power level in W for scrambling pulse
;zgoptns : Dfsat,Dfscuba,Dfsat_diff,Dscramble,DIz,Dw_gate