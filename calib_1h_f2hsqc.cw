;1H calibration by nutation during zz filter
; set pl18, applied at cnst18 (ppm)
;avance-version (05/10/28)
;HSQC
;2D H-1/X correlation via double inept transfer
;   using sensitivity improvement
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;constant time version
;using trim pulses in inept transfer
;using shaped pulses for inversion on f2 - channel
;
;(G.W. Vuister & A. Bax, J. Magn. Reson. 98, 428-435 (1992))
;(A.G. Palmer III, J. Cavanagh, P.E. Wright & M. Rance, J. Magn.
;   Reson. 93, 151-170 (1991))
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


;$OWNER=nmrsu
prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>



"p2=p1*2"
"p4=p3*2"
"p22=p21*2"
"d0=3u"
"d4=1s/(cnst2*4)"
"d11=30m"
"d12=20u"

"d18=0.1u"
"in18=inf1" ; nutation time increment

"d20=d23-p16-d16-p14*1.5-4u-d12"
"d31=d20"

"DELTA1=p16+d16+8u"
"DELTA2=d4-larger(p2,p8)/2-4u-p16-d16"
"DELTA3=d23-d0-p14/2-larger(p14,p22)-4u"
"DELTA4=d24-p19-d16"
"DELTA5=d4-larger(p2,p8)/2-4u-p16-d16"

"acqt0=0"

"spoff3=0"
"spoff5=bf2*(cnst21/1000000)-o2"
"spoff13=0"


1 ze 
  d11 pl12:f2 pl26:f3
2 d1 do:f2 do:f3


 4u pl1:f1 pl2:f2 pl3:f3

  ; first INEPT
  (p1 ph1)
  p16:gp5
  d16
  DELTA2 pl0:f2
  4u
  (center (p2 ph1) (p8:sp13 ph6):f2 )
  4u
  DELTA2 pl2:f2 UNBLKGRAD
  p16:gp5
  d16
  (p1 ph2)

; zz
  4u
  p16:gp6
  d16

; calibration
5u fq=cnst18(bf ppm):f1
5u pl18:f1
d18 cw:f1 ph1
5u do:f1
5u fq=0:f1 pl1:f1
p16:gp7  ; cleaning gradient
200u

; resume hsqc
  (p3 ph3):f2
  d0
  (center (p2 ph7) (p14:sp5 ph1):f2 (p22 ph1):f3 )  ; CO
  4u
  DELTA3 pl0:f2
  (p14:sp3 ph4):f2  ; Cali
  d20
  p16:gp1*EA*-1
  d16 pl0:f2
  (p14:sp5 ph1):f2   ; CO
  4u
  d12 pl2:f2

  (ralign (p1 ph1) (p3 ph4):f2 )
  p19:gp3
  d16
  DELTA4
  (center (p2 ph1) (p4 ph1):f2 )
  DELTA4
  p19:gp3
  d16
  (center (p1 ph2) (p3 ph5):f2 )
  p16:gp4
  d16
  DELTA4 pl0:f2
  (center (p2 ph1) (p8:sp13 ph1):f2 )
  DELTA4
  p16:gp4
  d16
  (p1 ph1)  ; final echo for gradient selection
  DELTA1
  (p2 ph1)
  4u
  p16:gp2
  d16 pl12:f2 pl26:f3
  4u BLKGRAD

  go=2 ph31 cpd2:f2 cpd3:f3
  d1 do:f2 do:f3 mc #0 to 2 
     F1QF(id18)
d31
exit 
  

ph1=0
ph2=1
ph3=0 2
ph4=0 0 2 2
ph5=1 1 3 3
ph6=0
ph7=0 0 2 2
ph31=0 2 2 0


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;pl26: f3 channel - power level for CPD/BB decoupling
;sp3 : f2 channel - shaped pulse 180 degree (on resonance)
;sp5 : f2 channel - shaped pulse 180 degree (off resonance)
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p14: f2 channel - 180 degree shaped pulse
;p16: homospoil/gradient pulse
;p19: 2nd homospoil/gradient pulse
;p22: f3 channel - 180 degree high power pulse
;p28: f1 channel - trim pulse                        [1 msec]
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d20 : = d23
;d23: d23 = T : 13.3 or 26.6 msec
;     2T (constant time period) = n/J(CC)
;d24: 1/(8J)XH for all multiplicities
;     1/(4J)XH for XH
;cnst2: = J(XH)
;cnst21: CO chemical shift (offset, in ppm)
;in0: 1/(2 * SW(X)) = DW(X)
;in20: 1/(2 * SW(X)) = DW(X) = in0
;nd0: 2
;NS: 4 * n
;DS: 32
;td1: number of experiments
;FnMODE: echo-antiecho
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    gp 1 : gp 2 : gp3 : gp4
;                         80 : 20.1 : 11  : -5   for C-13
;                         80 :  8.1 : 11  : -5   for N-15

;for z-only gradients:
;gpz1: 80%
;gpz2: 20.1% for C-13, 8.1% for N-15
;gpz3: 11%
;gpz4: -5%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SINE.100
;gpnam4: SINE.100



;$Id: hsqcctetgpsisp,v 1.4 2005/11/10 12:16:59 ber Exp $
