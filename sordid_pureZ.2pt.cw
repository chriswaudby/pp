;diffusion no D1, triple inept, 20-12-2010

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

prosol relations=<triple>

;define the number and intensity of gradients with respect to 100% [choose values between 0.1 and 0.9]
define list<gradient>  diff={0.15 1}
; incell:
;define list<gradient> diff={0.4 1}
;define list<gradient> diff={0.07 0.7 0.14 0.63 0.21 0.56 0.28 0.49 0.35 0.42}  ;10 points
;define list<gradient> diff={0.07 0.14 0.21 0.28 0.35 0.42 0.49 0.56 0.63 0.7}  ;10 points

define delay scanTime

;d20: 	total diffusion delay (first gradient first block to first gradient second block in next scan)
;cnst4: JNH
"d2 = 0.25/cnst4"
"d15 = 10u"
"d16 = 150u"
"d11 = 2m"
"p11=968245/(cnst21-cnst20)"  ; selective 90=sqrt(15)/(4*delta v)

"d17 = d25*0.5 - p21 - d16"
"d7 = d25*0.5 - p21 - d15"
"d4 = d2 - p30 - d16 - d15 - p11*0.637 - d7"
"d14 = d2 - p30 - d16 - d15 - p11*0.637 - d17"
"d8 = d2 - p16 - d16 - d15 - p21*0.637"

"d27=2*(5*d16+d4+2*p11+5*d15+p16*4+p1*4.83333+3u+p21*4+d8*2+d14)"
"d28=(p30*4+d16*4+d17*2+d7*2+p21*4+d15*5+d4+d21+d14+p11+p16+d11*2)"


;"d1 =d20-d11*4-d27-d28-aq" ; CW - ERROR in calculating big DELTA!
"d1 =d20-d11*6-d27-d28-aq-d15*10-d16*10-d2-d25-de-3m-p30*2-4u" ; calculate required recycle delay given d20

"d21=p11*0.637"
;"l4=l3/2"

"scanTime=d20-d15*2-d16*2-d2*4-p1*4.833-p11*0.726-p16*2-p21*4.726-5u" ; total time for single scan


;d20 : delay between first gradient first block and first gradient last block of next scan
;d25 : delay between the two gradients of the bipolar pairs

;"cnst20=2349"  ; H2O frequency in Hz (500MHz)
;"cnst21=4149"  ; center amide region in Hz (8.3ppm, 500MHz)

;p30 : duration of the encoding/decoding gradients


1 ze
  scanTime
  d11 igrad diff ;sync diff list ID with buffer ID
2 d11*2 
3 d11*2
4 d1 pl1:f1 pl3:f3
  d11
  d11
;************* first inept Nz(enc -1) -> Nz(enc -1), Hz -> 2HzNz(enc)
  d15 UNBLKGRAD
  d15
  p16:gp10  ; suppress residual Nx and Ny (if T2 comparable to DELTA)
  d16
  (p11:sp1 ph5):f1    	; selective 90 (square pulse on cnst21)
  d14 pl1:f1
  d15
  p30:gp6*diff                      ;gradient encoding
  d16
  d17
  ;(center (p12*0.83333:sp2 ph11 1u p12*3.16666667:sp2  ph12 1u p12*0.83333:sp2 ph11):f1 (p21*2 ph0):f3)  ;selective 180 on 1H
  (center (p1*0.83333 ph11 1u p1*3.16666667  ph12 1u p1*0.83333 ph11):f1 (p21*2 ph0):f3)  ;hard composite 180 on 1H
  d7
  d15
  p30:gp6*-1*diff                     ;gradient encoding
  d16
  d4
  (p11:sp1 ph3):f1 	; selective 90
  d15 pl1:f1
  p16:gp0
  d16
  (p1*0.83333 ph11 1u p1*3.16666667 ph12 1u p1*0.83333 ph11):f1
  3u
 
 ;************* second inept Nz(enc -1) -> 2HzNz(enc -1), 2HzNz(enc) -> Nz(enc) 
  (p21 ph7):f3
  d8 dgrad diff
  d15
  p16:gp1
  d16
  (center (p1*0.83333 ph11 1u p1*3.16666667 ph12 1u p1*0.83333 ph11):f1 (p21*2 ph0):f3)  ; hard composite 180 on 1H
  d15
  p16:gp1
  d16
  d8
  (p21 ph1):f3
  d15
  p16:gp2
  d16
   
 
 ;************* third inept 2HzNz(enc -1) -> Hz(enc -1)(dec -1), Nz(enc) -> Nz(enc) 
  (p11:sp1 ph4):f1    	;	selective 90
  d14 pl1:f1
  d15
  p30:gp6*diff                      ;gradient encoding
  d16
  d17
  (center (p1*0.83333 ph21 1u p1*3.16666667 ph22 1u p1*0.83333 ph21):f1 (p21*2 ph6):f3)  ;hard composite 180 on 1H
  d7
  d15
  p30:gp6*-1*diff              ;gradient encoding
  d16
  d4*0.5 igrad diff
  d4*0.5 igrad diff 
  d21 BLKGRAD
 ; aqcuisition
  goscnp ph31
  d11 st  ; increment buffer
  d11 ipp4 ipp7 ipp21 ipp22 ipp31
  lo to 2 times nbl  ; inner loop (through diff list)
  d11 
  d11 
  lo to 3 times ns  ; outer loop
  d11 wr #0  ; write data to disk
  d11 zd
  lo to 4 times td0  ; repeat whole acquisition block
exit
 
 ph0=0
 ph1=1
 ph2=2
 ph3=3
 ph4=2 2 2 2 2 2 2 2
     3 3 3 3 3 3 3 3
     0 0 0 0 0 0 0 0
     1 1 1 1 1 1 1 1
 ph5=0
 ph6=0 
 ph7=0 2 0 0
 ph11=(360) 34  ;304  ;ph11 and ph12 for a y pulse (test whether should not be negative)
 ph12=(360) 144  ;54
 ph21=(360) 34 34 34 34 124 124 124 124 214 214 214 214 304 304 304 304
 ph22=(360) 144 144 144 144 234 234 234 234 324 324 324 324 54 54 54 54
 ph31=0 2 2 0 2 0 0 2
      3 1 1 3 1 3 3 1
      2 0 0 2 0 2 2 0
      1 3 3 1 3 1 1 3
 

;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;sp1: f1 channel - shaped pulse 90 degree (SQUARE)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p11: f1 channel -  90 degree selective excitation
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse
;p30: encoding gradient (delta/2)
;d20: diffusion delay (big DELTA)
;cnst4: = J(YH)
;cnst20:= H2O frequency in Hz
;cnst21:= center of amide frequency in Hz
;NS: 1 * n
;DS: >= 16
;td1: length of gradient list
;FnMODE: QF

;for z-only gradients:
;gpz0: 13.17%
;gpz1: 11.73%
;gpz2: 17.13%
;gpz6: Gmax%
;gpz10:15%

;use gradient files:   
;gpnam0: SINE.10
;gpnam1: SINE.10
;gpnam10 SINE.10
;gpnam2: SINE.10
;gpnam6: SMSQ10.100
