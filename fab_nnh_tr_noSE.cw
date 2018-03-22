;zg
;avance-version
;1D sequence
;for cryoprobes introduced p20 pulse for h2o purge 2ms
;be careful when using in conv. probe use 3D gradients

#include <Avance.incl>
#include <Grad.incl>

;#define ONE_D


"d16=200u"
"d15=5u"
"d2 = 2.7m"
"d4 = d2-p16-d16-d15"
"d5 = d2-p16-d16-d15-p11"
"d30 = 30m"
"d0 = 3u"
"d10 = p1*2+d0*2"
"d22 = d2+p1+d0"
"d23 = d2-p1-d0-3u"
"d24 = d20*0.5-d15-p16-d16-p2-p2*0.637"


"l5=1"
"in0=inf1/2"

1 ze
2 d30 do:f3
3 d1 pl1:f1 pl3:f3 
  10u
  (p2 ph0):f3
  d15 UNBLKGRAD
  p16:gp0
  1m 
  (p1 ph0):f1 
  d15 
  p16:gp1
  d16  
  d4
  (center (p1*2 ph2):f1 (p2*2 ph0):f3)
  d4
  d15
  p16:gp1
  d16
  (p1 ph4):f1  
  d15
  p20:gp2
  1m

if "l5 %4 == 1" goto 20
if "l5 %4 == 2" goto 20

  (p2 ph10):f3
  d2
  (center (p1*2 ph0):f1 (p2*2 ph0):f3)
  d2
  (p2 ph3):f3
  d15
  p16:gp5
  d16
  (p1 ph0):f1
  1u

 goto 30

20  d15
  p16:gp5
  d16
  (p1 ph0):f1
  3u
  (p1 ph10):f1
  1u
  
;******* relaxation *******
30 3u 
  (p2 ph11):f3
  d24
  d15 
  p16:gp6
  d16 
  (p2*2 ph0):f3
  d15
  p16:gp6
  d16
  d24
  (p2 ph12):f3
  3u

;******** end relaxation **********

if "l5 %4 == 1" goto 40
if "l5 %4 == 3" goto 40
  
  1u
  (p1 ph0):f1
  d15
  p16:gp7
  d16
  (p2 ph5):f3
  d22
  (p2*2 ph24):f3
  d0
  (p1*2 ph1):f1 
  d23
  d0
  (p2 ph3):f3
  3u
  goto 50

40 1u
  (p1 ph0):f1
  3u
  (p1 ph14):f1
  d15
  p16:gp7
  d16
  (p2 ph5):f3
  d0
  (p1*2 ph0):f1 
  d0
  (p2*2 ph0):f3
  d10 
  (p2 ph0):f3

50  d15 
   p16:gp3
   d16
   (p1 ph0):f1  
   d15
   p16:gp4
   d16
   d5
   3u pl0:f1
   3u
  (p11:sp1 ph23:r):f1
   3u
   3u pl1:f1
   3u
   (center (p1*2 ph1):f1 (p2*2 ph0):f3)
  3u
  3u pl0:f1
  3u
  (p11:sp1 ph23:r):f1
  6u
  d15 
  p16:gp4
  d16 BLKGRAD
  d5 pl16:f3
  go=2 ph31 cpds3:f3
  d30 do:f3 mc #0 to 2
     F1I(iu5, 4)
     F1PH(ru5 & ip5, id0)
  d30 do:f3
exit

ph0=0 
ph1=1
ph2=2
ph3=3
ph23=3
ph4={1}*8 {3}*8
ph5={0}
ph10 = {0}*4 {2}*4
ph11 = {0}*1 {2}*1
ph12 = {0}*16 {2}*16
ph14 = {0}*2 {2}*2
ph15 = {2}*2 {0}*2
ph20 =  {2}*4 {0}*4
ph24 = {0}*2 {1}*2


ph31=0 2 2 0 2 0 0 2
     2 0 0 2 0 2 2 0
     2 0 0 2 0 2 2 0
     0 2 2 0 2 0 0 2
     


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  high power pulse
;d1 : relaxation delay; 1-5 * T1
