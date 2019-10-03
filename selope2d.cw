;SELOPE_2D.js
;2D 1H-1H correlation experiment using selective excitation and J-tranfer 
;J. Schlagnitweit, E. Steiner, H. Karlsson, K. Petzold
;2017


prosol relations=<triple>

#include <Avance.incl> 
#include <Grad.incl> 
#include <Delay.incl>

"p2=p1*2"
"d12=20u"

;"spoff4=bf1*(cnst21/1000000)-o1" 
"spoff4=0"
"cnst28=0"

"TAU=de+p1*2/3.1416+50u"
"acqt0=0"
baseopt_echo

"d0=3u"
"in0=inf1"



1 ze 
2 30m

  d1 fq=cnst28:f1  ; on-resonance on signals of interest
  50u ;UNBLKGRAD

;  1u
;  p16:gp5
;  d16
; 1u pl0:f1
; (p13:sp4 ph23):f1 ; 1u pl1:f1
; p1 ph10
; 1u
; p16:gp3
; d16 pl0:f1

(p13:sp4 ph9):f1

d0

d5 pl1:f1  ; INEPT
p2 ph0
d5
p1 ph12 ; y
d5
p2 ph0 
d5

; excitation sculpting
1u fq=cnst29(bf ppm):f1 ; switch to H2O
50u UNBLKGRAD
p16:gp1
d16 pl0:f1
(p12:sp1 ph2:r):f1
4u
d12 pl1:f1

p2 ph3

4u
p16:gp1
d16
TAU
p16:gp2
d16 pl0:f1
(p12:sp1 ph4:r):f1
4u
d12 pl1:f1

p2 ph5

4u fq=cnst28:f1  ; back to o1
p16:gp2
d16 BLKGRAD

go=2 ph31
;30m mc #0 to 2 F1PH(calph(ph23, +90) & calph(ph10, +90) &calph(ph11, +90) &calph(ph21, +90) &calph(ph22, +90) &calph(ph9, +90), caldel(d0, +in0))
  30m mc #0 to 2 F1PH(calph(ph9, +90), caldel(d0, +in0))
  4u BLKGRAD

exit


ph0=0
ph1=0 0 0 0 2 2 2 2
;ph23=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph2=0 1
ph3=2 3
ph4=0 0 1 1
ph5=2 2 3 3
ph9=0 0 0 0 2 2 2 2 
ph12=3
ph31=0 2 2 0 2 0 0 2 ;2 0 0 2 0 2 2 0

;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;sp1 : f1 channel - shaped pulse 180 degree (ex. sculpting)
;sp4 : f1 channel - 90 sel pulse power level
;p1 : f1 channel - 90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p12: f1 channel - 180 degree shaped pulse for ex. sculpting (Squa100.1000)
;p13: f1 channel - 90 selective pulse (e.g. H8/H6 region) (Eburp2.1000)
;p16: homospoil/gradient pulse
;d1 : relaxation delay; 1-5 * T1
;d12: delay for power switching
;d16: delay for homospoil/gradient recovery 
;d5: INEPT delay ~25ms, use 1/4J
;ns: 8 * n, total number of scans: NS * TD0 
;ds: 4
;o1p: on-resonance for selective excitation
;cnst29: chem. shift [ppm] of water signal

;for z-only gradients: 
;gpz1: 31%
;gpz2: 11%
;gpz3: 49%
