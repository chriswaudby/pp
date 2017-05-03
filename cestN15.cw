;15N-CEST experiment based on Lakomek sensitivity-enhanced 15N-T1 measurement
;Chris Waudby, June 2016
;


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

;list of CEST saturation frequencies
;first line in file 'Hz sfo'
;second line in file should be zero, indicating the reference plane
;give values in Hz relative to sfo3
define list<frequency> N15sat = <$FQ3LIST>

;38 points in steps of 45 Hz:
;define list<frequency> N15sat = {0  -810.0 -765.0 -720.0 -675.0 -630.0 -585.0 -540.0 -495.0 -450.0 -405.0 -360.0 -315.0 -270.0 -225.0 -180.0 -135.0 -90.0 -45.0 0.0 45.0 90.0 135.0 180.0 225.0 270.0 315.0 360.0 405.0 450.0 495.0 540.0 585.0 630.0 675.0 720.0 765.0 810.0}
;define list<frequency> N15sat = {0 -1000 0}
;define list<frequency> N15sat = {0}


"in0=inf2*0.5"

"p2=p1*2"
"p22=p21*2"

"d0=10u"
"d11=30m"
"d26=p21-p1"
"DELTA=2.65m"
"DELTA1=2.65m"
"DELTA2=2.65m-200u"
"DELTA3=2.65m-200u"
"DELTA4=p24-p19*0.6366"
"DELTA5=2.65m-p25-200u"
"DELTA7=2.65m"

"l2=1" ; loop counter for saturation list
; used to identify first point as the reference spectrum,
; and to activate temperature compensation for this point

"DELTA6=30u+p2" ; for zero phase correction

;"cnst18=-800"
;"cnst19=8.2"


aqseq 312


1 ze
  1m
2 2m do:f3
  1m LOCKH_OFF
  d11
3 3m
4 1m
5 1m do:f3
  10u UNBLKGRAD

;---------temperature compensation and d1 recovery delay---------

; apply temperature compensation immediately after 1st point (reference with no spinlock)
if "l2==2" ; i.e. if following the first FID
{
    10u LOCKH_ON
    10u fq=cnst19(bf ppm):f1 ; put 1H on amides
    10u pl8:f1
    d18 cpds1:f1
    10u do:f1
    10u fq=0:f1 ; put 1H back on water
    10u LOCKH_OFF
}
; purge water before recycle delay
10u pl1:f1
(p1 ph0):f1
5u
p16:gp8
200u
1m BLKGRAD

d1

1m UNBLKGRAD

;------- kill steady state 15N ------------
10u pl3:f3
(p21 ph0):f3
5u
p16:gp6
200u

;------- first INEPT Hz-> 2HxNz -------------------------
(p1 ph0):f1
5u
DELTA gron0 ; soft gradient to prevent radiation damping
5u  groff
(center (p2 ph0):f1 (p22 ph0):f3)
5u
DELTA gron0
5u groff

;------- rephase  2HxNz to Nz-----------------------------
(center (p1 ph5):f1  (p21 ph0):f3)
5u
DELTA1 gron1 ; soft gradient to prevent radiation damping
5u groff
(center (p2 ph0):f1 (p22 ph0):f3)
5u
DELTA1 gron1
5u groff
(center (p1 ph2):f1 (p21 ph6):f3) ; phase-cycle Nz, -Nz for Freeman-Hill decay
                                  ; purge pulse to kill any residual HzNz
;--------------------------------------------------------
5u
p16:gp6  ; cleaning gradient
200u

;------15N T1 relaxation period--------------------------
5u pl8:f1 pl18:f3

if "l2==1" goto 77

5u LOCKH_ON
5u N15sat:f3
5u fq=cnst19(bf ppm):f1
d18 cpds1:f1 cw:f3 ph0
5u do:f1 do:f3
5u LOCKH_OFF

77 5u pl1:f1 pl3:f3
   5u fq=0:f1
   5u fq=0:f3

5u
p16:gp7  ; cleaning gradient
200u

;---------- Echo-Antiecho encoding------------------------
(p21 ph7):f3
DELTA6 ;compensation for d0 15N evolution
DELTA5
190u
p25:gp5*EA ; 1ms
10u
(center (p2 ph0):f1 (p22 ph7):f3)
10u
p25:gp5*EA*-1 ; 1ms
190u
DELTA5

; t1 evolution ---------------------------------------------
d0*0.5 gron1
5u groff
d0*0.5 gron1*-1
5u groff
(p2 ph0):f1
d0
;--------Rance-Kay transfer back ----------------------------

(center (p1 ph0):f1 (p21 ph8):f3)
DELTA2 gron2
200u groff
(center(p2 ph0):f1 (p22 ph0):f3)
DELTA2 gron2
200u groff

; ---second INEPT ----------------------------------------
(center (p1 ph1):f1 (p21 ph1):f3)    ;DOUBLE 90
DELTA3 gron3
200u groff
(center(p2 ph0):f1 (p22 ph0):f3)
DELTA3 gron3
200u groff

d26*0.5 ;d26 = p21 - p1
(p19 ph0:r):f1 ; fine tune p19 for water suppression
245u
DELTA4
(p2 ph0):f1
5u
p24:gp4 ; 201u
200u

999  10u
10u pl31:f3
20u BLKGRAD

go=2 ph31 cpds3:f3
500u do:f3
d11 wr #0 if #0 zd

; loop over CEST frequencies
1m N15sat.inc
1m iu2
lo to 3 times l6

; Echo/Anti-echo
1m ip8*2
1m igrad EA
1m ru2
1m N15sat.res
lo to 4 times 2

; t1 evolution
1m id0
lo to 5 times l3

1m do:f3
1m BLKGRAD
exit

ph0= 0
ph1= 1          ;check right phase for Boltzmann !!!!!
ph2= 2
ph5= 1 1 1 1
ph6= 1 1 3 3
ph7= 1 3
ph8= 0
ph31=1 3 3 1


;-------------NOTES----------------------

;o1p = 4.7 ppm
;o3p=119 ppm

;NS=4*n
;in0=inf/2
;SW=1/(2*in0)
;echo-antiecho in N15 (process as Complex in NmrDraw before splitting the spectra)

;d18: relaxation time

;l6: number of CEST frequencies (including reference) (td1)
;l3: number of complex points (td2 / 2)
;cnst19: chemical shift at center of amide region (ppm)


; 1H pulses

;p1:  90 deg hard 1H pulse @pl1
;p19 last 90 deg hard 1H pulse @pl1, can be adjusted for improving water-supression
;pl1: 1H 90 deg
;pl0: 120 dB
;pcpd1: 1H 90 for CPD (4.5 kHz at 900 MHz = 55.5 us)
;pl8: 1H decoupling power (4.5 kHz at 900 MHz = 55.5 us)


;15N pulses
;p21 : 90 deg hard 15N pulse @pl3
;pl3 :15N 90 deg
;pl18: 15N spin-lock power
;CPDPRG3: garp (aq 15N decoupling)
;pcpd3: 15N decoupling (200u @pl31)
;pl31: 15N decoupling power



; gradients
;p16: 1000u
;p17: 200u
;p24: 201u Echo/Anti-echo decoding gradient
;p25: 1000u Echo/Anti-echo half-encoding gradient

;for z-only gradients
;gpz10: -1.5%
;gpz2: -5%
;gpz3: -0.5%
;gpz4: 40%
;gpz5: 40%
;gpz0: 3%
;gpz1: 2%
;gpz6: 30%
;gpz7: -50%
;gpz8: 77%

;gpnam4 SINE.10
;gpnam5 SINE.10
;gpnam6 SINE.50
;gpnam7 SINE.10
