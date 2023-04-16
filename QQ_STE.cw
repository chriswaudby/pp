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

; pulses
"p2=p1*2"
"p4=p3*2"
"p31=p30*0.5"

; delays
"in0=inf1/2"
"d11=30m"

define list<gradient> diff=<Difframp>

define delay taua
"taua=d5"   ; < 1/(4JCH) 1.8 ms

define delay taub
"taub=d6"   ; = 1/(4JCH) 2.0 ms

define delay T_diff
"T_diff = d20"

"acqt0=-p1*2.0/PI - 2u"  ; select 'DIGIMOD = baseopt' to execute

1 ze

2 d11 do:f2


; presat during D1
20u fq=cnst1:f1         ; jump from o1 to water
4u pl10:f1              ; power(tsatpwr) for presaturation
d1 cw:f1 ph26           ; Hcw(d1)x
4u do:f1                ; cw off
2u pl1:f1               ; power(tpwr)


2u pl1:f1
2u pl2:f2

20u UNBLKGRAD   ; dly 20u, unblank gradients and lock on hold
10u fq=0:f1     ; jump back to o1

(p3 ph26):f2   ; C90x - To destroy 13C Boltzman magnetization
  2u
  (p16:gp0)     ; gradient 0
  d16


; This is the real start

(p1 ph1):f1                    ; H90ph1

2u
p16:gp1                         ; gradient 1
d16
"DELTA = taua - 2.0u - p16 - d16 - p3 - p1*4/PI - 0.2u"
DELTA                                    ; delay 1/4JCH

(center (p2 ph1):f1 (p4 ph26):f2) ; H180x,C180x

"DELTA = taua - 2.0u - p16 - d16 - p3 - 2u"
DELTA                                    ; delay 1/4JCH
2u
p16:gp1           ; gradient 1
d16
2u pl2:f2

(p3 ph3):f2      ; C90ph3

2u
p16:gp2          ; gradient 2
d16
"DELTA = taub - 2.0u - p16 - d16"
DELTA             ; delay 1/4JCH

(center (p2 ph1):f1 (p3 ph27 p4 ph26 p3 ph27):f2) ; H180x,C180x

"DELTA = taub - 2.0u - p16 - d16"
DELTA             ; delay 1/4JCH
2u
p16:gp2          ; gradient 2
d16

(p1 ph1):f1      ; H90ph1 - FIRST DIFFUSION ENCODING PERIOD

2u
p31:gp6*diff      ; gradient 6 diffusion encoding gradient
d17

(p3 ph27 p4 ph26 p3 ph27):f2   ; C180x

2u
p31:gp6*diff      ; gradient 6 diffusion encoding gradient
d17

(p1 ph27 p2 ph26 p1 ph27):f1

2u
p31:gp6*-1.0*diff ; -gradient 6 diffusion encoding gradient
d17

(p3 ph27 p4 ph26 p3 ph27):f2   ; C180x

2u
p31:gp6*-1.0*diff ; -gradient 6 diffusion encoding gradient
d17

; transfer back to longitudinal for diffusion delay
(p1 ph26):f1          ; H90x

2u
p16:gp3               ; gradient 3
d16
"DELTA = taub - 2.0u - p16 - d16"
DELTA                 ; delay 1/4JCH

(center(p1 ph27 p2 ph26 p1 ph27):f1 (p3 ph27 p4 ph26 p3 ph27):f2) ; H180x,C180x

"DELTA = taub - 2.0u - p16 - d16"
DELTA                 ; delay 1/4JCH
2u
p16:gp3               ; gradient 3
d16

(lalign (p1 ph27):f1 (p3 ph26):f2)         ; H90y, C90x


; magnetization diffuses as 2IzCz
"DELTA = T_diff - taub*4 - 2u - p16 - d16"
DELTA
2u
p16:gp4                 ; gradient 4
d16


; transfer back to 4Q for decoding

(ralign (p1 ph27):f1 (p3 ph26):f2)          ; H90y, C90X

2u
p16:gp5               ; gradient 5
d16
"DELTA = taub - 2.0u - p16 - d16"
DELTA                 ; delay 1/4JCH

(center(p2 ph26):f1 (p3 ph27 p4 ph26 p3 ph27):f2) ; H180x,C180x

"DELTA = taub - 2.0u - p16 - d16"
DELTA                 ; delay 1/4JCH
2u
p16:gp5               ; gradient 5
d16

(p1 ph26):f1          ; H90ph6

2u
p31:gp6*diff      ; gradient 6 diffusion encoding gradient
d17

(p3 ph27 p4 ph26 p3 ph27):f2   ; C180x

2u
p31:gp6*diff      ; gradient 6 diffusion encoding gradient
d17

(p1 ph27 p2 ph26 p1 ph27):f1    ; H180 composite

2u
p31:gp6*-1*diff      ; gradient 6 diffusion encoding gradient
d17

(p3 ph27 p4 ph26 p3 ph27):f2   ; C180x

2u
p31:gp6*-1*diff      ; gradient 6 diffusion encoding gradient
d17

(p1 ph2):f1          ; H90ph2

2u
p16:gp7               ; gradient 7
d16

"DELTA = taub - 2.0u - p16 - d16"
DELTA

(center (p1 ph27 p2 ph26 p1 ph27):f1 (p3 ph27 p4 ph26 p3 ph27):f2) ; H180x,C180x

"DELTA = taub - 2.0u - p16 - d16"
DELTA

2u
p16:gp7                ; gradient 7
d16

(p3 ph26):f2        ; C90x

2u
p16:gp8             ; gradient 10
d16

"DELTA = taua - 2.0u - p16 - d16 - p3 + p1*2/PI"
DELTA   ; delay 1/4JCH

(center(p1*2 ph26):f1 (p3*2 ph26):f2) ; H180x,C180x

2u pl2:f2
"DELTA = taua - 2u - 2u - p16 - d16 - p3 - 4u" ; - p3 - p3 - p1*2/PI"
DELTA                                    ; delay 1/4JCH

2u
p16:gp8
d16

4u BLKGRAD
;(p3 zero):f2
;(p3 ph4):f2

;(p1 one):f1


2u pl21:f2

go=2 ph31 cpds2:f2

d11 do:f2 mc #0 to 2 F0(zd) ; write fid to disk
  F1QF(calgrad(diff))

HaltAcqu, 1m
exit


ph1= (6) 0 1 2 3 4 5  ; 0 60 120 180 240 300
ph2= (6) 0 0 0 0 0 0
         1 1 1 1 1 1
         2 2 2 2 2 2
         3 3 3 3 3 3
         4 4 4 4 4 4
         5 5 5 5 5 5
ph3= 0
ph26=0
ph27=1
ph28=2
ph29=3
ph31=(6) 0 3 0 3 0 3
         4 1 4 1 4 1
         2 5 2 5 2 5
         0 3 0 3 0 3
         4 1 4 1 4 1
         2 5 2 5 2 5


;d1 : repetition delay
;d5 : taua ~1/4JCH
;d6 : taub = 1/4JCH 2 ms
;d11 : delay for disk i/o, 30ms
;d16 : gradient recovery delay, 200us
;d19 : delay for binomial water suppression = |1/(4*dis)| dis = Hz between o1 and water
;      d19 = (1/(4*|cnst1|)), cnst1 = distance between methyl and water (in Hz)
;d20: T_diff total diffusion time
;pl1 : tpwr - power level for p1
;pl2 : dhpwr - power level for hard 13C pulse p3
;pl21 : dpwr - power level for  13C decoupling cpd2
;p1 : p1
;p3 : p3
;cpd2 : 13C decoupling according to program defined by cpdprg2
;pcpd2:  13C 90 degree pulse at pl21 for cpd2
;spnam23 : file name for chirp pulse during inepts
;p16 : duration of diffusion encoding gradient pulse
;GPZ3 : set to 100 % for gradient encoding
;cns31 6 water - 01 (Hz)
;cnst12 : power level in W for scrambling pulse
;zgoptns : Dfsat,Dfscuba,Dfsat_diff,Dscramble,DIz,Dw_3136