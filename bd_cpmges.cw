;bd_cpmges
;avance-version (00/10/05)
;excitation sculpted 1D with CPMG filter
; Modified by Pedro M. Aguiar
; 2015 10 08
;Modified chris waudby, Mar 2020 - correcting CPMG period, using baseopt


; set d21 so that d20 > 0
; d21 should be << 1/J and > 40*p2 [0.5 to 2.0 ms]
; use loop counter l1 to set duration of cpmg filter period 


prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>


"p2=p1*2"
"d12=20u"
"p5=200u"

"d20=((d21/2) - d16 - p5/1e6 - p1/1e6)"

define delay cpmg_time

"cpmg_time = l1*d21"

"acqt0=-0.6366*p1"

cpmg_time



1 ze
2 d12 BLKGRAD
  d1 pl1:f1 UNBLKGRAD
  p1 ph1

3 p5:gp3
  d16
  d20
  (p2 ph6):f1
  p5:gp3
  d16
  d20
  lo to 3 times l1
  
  50u 
  p16:gp1
  d16 pl0:f1
  (p12:sp1 ph2:r):f1
  4u
  d12 pl1:f1

  p2 ph3

  4u
  p16:gp1
  d16 
  50u
  p16:gp2
  d16 pl0:f1
  (p12:sp1 ph4:r):f1
  4u
  d12 pl1:f1

  p2 ph5

  4u
  p16:gp2
  d16

  go=2 ph31
  wr #0
  4u BLKGRAD
exit


ph1=0
ph2=0 1
ph3=2 3
ph4=0 0 1 1
ph5=2 2 3 3
ph6=1
ph31=0 2 2 0 


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;sp1 : f1 channel - shaped pulse 180 degree
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p12: f1 channel - 180 degree shaped pulse (Squa100.1000)   [2 msec]
;p5 : gradient pulses in CPMG loop (fixed to 200us)
;p16: homospoil/gradient pulse
;d1 : relaxation delay; 1-5 * T1
;d12: delay for power switching                             [20 usec]
;d16: delay for homospoil/gradient recovery
;d20: set automatically using d21 
;d21: half duration of each CPMG loop [0.5-2.0 ms]
;NS: 8 * n
;DS: 4
 


;use gradient ratio:    gp 1 : gp 2
;                         31 :   11

;for z-only gradients:
;gpz1: 31%
;gpz2: 11%
;gpz3: 8%

;use gradient files:
;gpnam1: SINE.100
;gpnam2: SINE.100



;$Id: zgesgp,v 1.3 2000/10/06 09:09:37 ber Exp $
