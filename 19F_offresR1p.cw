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
| adiabatic_flg: use adiabatic passage to SL angle
|
| Pseudo-3D
| Jan Overbeck
| 2020
|
*/

/*--------------------------------
;     Parameters to set
; -------------------------------*/
;p11: adiabatic pulse length, 4ms
;spnam4: adiabatic ramp = tanhtan90
;sp4: adiabatic ramp power, = pl25
;pl25: spin lock power, = sp4
;p32: spin lock lenght Tex
;VPLIST: list of spin lock lengths
;FQ1LIST: list of spin lock offsets !bf hz!

#include <Avance.incl>
#include <Grad.incl>

define list<pulse> plength = <$VPLIST>
define list<frequency> fqlist = <$FQ1LIST>

"p2=p1*2"
"d11=30m"
"l2=0"
"l3=0"
"cnst28=fqlist"

#ifndef adiabatic_flg
"p3 = p1*pow(10,(10*log10(plw1) – 10*log10(plw25))/20)"       ;90 degree SL pulse
"p6 = ((cnst28)/((1/(p3*4))))"               ; spin lock offset / spin lock power
"p7 = atan(p6)"                          ; arc tan from this ratio = angle in rad
"p8 = p7*360/(2*PI)"                                            ; angle in degree
"p4 = p1*(1-p8/90)"                                            ; new pulse length
#endif

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
  "p7 = atan(p6)"                       ; arc tan from this ratio = angle in rad
  "p8 = p7*360/(2*PI)"                                         ; angle in degree
  "p4 = p1*(1-p8/90)"                                         ; new pulse length
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

/*    cycle above&below plane     */
;  1u fq=cnst29(bf ppm):f1
;  1u pl1:f1
;  p1 ph7
;  p1 ph8
  if "p32 == 0.0"
  {
    30m
    #ifdef adiabatic_flg
      1u fq=fqlist:f1
      1u pl25:f1
      (p11:sp4(currentpower) ph1):f1
    #else
      1u fq=cnst29(bf ppm):f1
      1u pl1:f1
      p4 ph4
    #endif
  }
  else {
    30m
    #ifdef adiabatic_flg
      1u fq=fqlist:f1
      1u pl25:f1
      (p11:sp4(currentpower) ph6):f1
    #else
      1u fq=cnst29(bf ppm):f1
      1u pl1:f1
      p4 ph4
    #endif
    1u fq=fqlist:f1
    1u pl25:f1
    (p32 ph6):f1                                    ; <-- this is the Spin Lock
  }
; ----------------------------------

/* ---------------------------------
;  transfer back to z
; --------------------------------*/
#ifdef adiabatic_flg
    1u fq=fqlist:f1
    1u pl25:f1
    (p11:sp5(currentpower) ph6):f1
#else
    1u fq=cnst29(bf ppm):f1
    1u pl1:f1
    p4 ph5
#endif
; ----------------------------------

/*  cycle above&below plane BACK  */
;  1u fq=cnst29(bf ppm):f1
;  1u pl1:f1
;  p1 ph7
;  p1 ph8

/* ---------------------------------
;     anti-ringing
; --------------------------------*/
  1u pl1:f1
  1u fq=cnst29(bf ppm):f1
  p1 ph1

    d13
  p1 ph2
  d13
  p1 ph3
; ----------------------------------

;  4u BLKGRAD
  go=2 ph31
  30m mc #0 to 2
  F1QF(calclc(l2,1))
  F2QF(calclist(fqlist,1))
exit

HaltAcqu, 1m
exit

ph1=0
ph2=2 2 0 0
ph3=0 0 0 0 2 2 2 2 1 1 1 1 3 3 3 3
ph4=1 3
ph5=3 1
ph6=0 2
ph7=0
ph8=0 2
ph31=0 0 2 2 2 2 0 0 1 1 3 3 3 3 1 1

;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O    [30 msec]
;ns: 16 * n
;ds: > 128
