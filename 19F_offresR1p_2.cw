/*
|
| Off-resonance 19F R1rho as pseudo-3D
|
| with different SL lenghts read in via VPLIST
| and different SL offsets read in via FQ1LIST
|
| VPLIST: spin lock lengths
| FQ1LIST: offsets in Hz
|
| Pseudo-3D
| Chris Waudby
| 2025
|
*/

/*--------------------------------
;     Parameters to set
; -------------------------------*/
;p44: f1 channel - 180 degree shaped pulse
;sp30: f1 channel - shaped pulse 180 degree (Bip720,50,20.1)
;pl25: spin lock power
;p32: spin lock lenght Tex
;VPLIST: list of spin lock lengths
;FQ1LIST: list of spin lock offsets !bf hz!

#include <Avance.incl>
#include <Grad.incl>

define list<pulse> plength = <$VPLIST>
define list<frequency> fqlist = <$FQ1LIST>

"p2=p1*2"
"d11=30m"
"l1=0"  ; z/-z and theta/theta-bar
"l2=0"  ; vplist
"l3=0"
"cnst28=fqlist"

"p29=10*log10(plw25)"
;"p3 = p1*pow(10,(10*log10(plw1) - 10*log10(plw25))/20)"       ;90 degree SL pulse
"p3 = p1*pow(10,(10*log10(plw1) - p29)/20)"       ;90 degree SL pulse
;"p3=9u"
"p6 = ((cnst28)/((1/(p3*4))))"            ; spin lock offset / spin lock power
"p7 = atan(p6)*180/PI"                    ; arc tan from this ratio = theta in deg
"p4 = p1*(1-p7/90)"                                          ; theta pulse length
"p5 = p1*(1+p7/90)"                                          ; 180-theta pulse length

aqseq 312

1 ze

/*--------------------------------
; calculate hard pulse for offset
; dependent tip angle theta
; -------------------------------*/
2 30m
  "p30 = plength.max"
  "p32=plength[l2]"
  "p31=p30-p32"
  "cnst28=fqlist"
  "p6 = ((cnst28)/((1/(p3*4))))"            ; spin lock offset / spin lock power
  "p7 = atan(p6)*180/PI"                    ; arc tan from this ratio = theta in deg
  "p4 = p1*(1-p7/90)"                                          ; theta pulse length
  "p5 = p1*(1+p7/90)"                                          ; 180-theta pulse length
; --------------------------------

/* ---------------------------------
;     heating compensation
; --------------------------------*/
if "p31 > 0.0"
  {
  1u fq=100(bf ppm):f1
  1u pl25:f1
  (p31 ph1):f1
  ;print "heating compensation on"
  }
; ----------------------------------

   d1
;50u UNBLKGRAD

/* ---------------------------------
;  transfer to theta and SL
; --------------------------------*/

1u fq=0:f1

if "l1 % 4 == 0"
{
; +z, theta
(p44:sp30 ph1):f1
1u
(p44:sp30 ph10):f1
1u pl1:f1
p4 ph2          ; theta(y)
1u fq=fqlist:f1
1u pl25:f1
(p32 ph1):f1    ; SL(x)
1u fq=0:f1
1u pl1:f1
p4 ph4          ; theta(-y)
}
if "l1 % 4 == 1"
{
; -z, theta
(p44:sp30 ph1):f1
1u pl1:f1
p4 ph4          ; theta(-y)
1u fq=fqlist:f1
1u pl25:f1
(p32 ph3):f1    ; SL(-x)
1u fq=0:f1
1u pl1:f1
p4 ph2          ; theta(y)
1u
(p44:sp30 ph10):f1
}
if "l1 % 4 == 2"
{
; +z, theta-bar
(p44:sp30 ph1):f1
1u
(p44:sp30 ph10):f1
1u pl1:f1
p5 ph2          ; thetabar(y)
1u fq=fqlist:f1
1u pl25:f1
(p32 ph3):f1    ; SL(-x)
1u fq=0:f1
1u pl1:f1
p5 ph4          ; thetabar(-y)
}
if "l1 % 4 == 3"
{
; -z, theta-bar
(p44:sp30 ph1):f1
1u pl1:f1
p5 ph4          ; thetabar(-y)
1u fq=fqlist:f1
1u pl25:f1
(p32 ph1):f1    ; SL(x)
1u fq=0:f1
1u pl1:f1
p5 ph2          ; thetabar(y)
1u
(p44:sp30 ph10):f1
}

/* ---------------------------------
;     anti-ringing
; --------------------------------*/
  1u pl1:f1
  p1 ph11
  d13
  p1 ph12
  d13
  p1 ph13
; ----------------------------------
  go=2 ph31
  30m mc #0 to 2
   F1I(iu1,4)
   F1QF(calclc(l2,1))
   F2QF(calclist(fqlist,1))
exit


ph1=0
ph2=1
ph3=2
ph4=3
ph10=2
ph11=0
ph12=2 0
ph13=0 0 2 2 1 1 3 3
ph31=0 2 2 0 1 3 3 1

;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O    [30 msec]
;ns: 16 * n
;ds: > 128