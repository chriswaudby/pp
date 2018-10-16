;Explicit coding of acquisition with prior setting of receiver phase
;Replaced trim pulse with ZZ period
;
;hsqcetgp
;avance-version (07/04/04)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;using trim pulses in inept transfer
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <De.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p4=p3*2"
"d4=1s/(cnst2*4)"
"d11=30m"
"d13=4u"


"d0=3u"
"in0=inf1/2"

"DELTA1=d4-p16-d16-4u-larger(p2,p8)*0.5"
"DELTA2=2*d0+p2"
"DELTA3=d16+0.6366*p1"


1 ze
  d11 pl1:f1 pl12:f2
2 d11 do:f2

  d1
  50u UNBLKGRAD
  4u pl1:f1 pl2:f2
 
   ; begin HSQC
  (p1 ph1)
  4u
  p16:gp1
  d16
  DELTA1
  (center (p2 ph2) (p8:sp13 ph1):f2 )
  DELTA1 
  p16:gp1
  d16 pl2:f2
  4u
  (p1 ph2)

  ; zz filter
  4u
  p16:gp2
  d16

  ; t1 evolution
  (p3 ph11):f2
  d0
  (p2 ph13)
  d0
  p16:gp3*EA
  d16
  (p4 ph4):f2
  p16:gp3*EA*-1
  d16
  DELTA2
  (p3 ph12):f2

  ; zz filter
  4u
  p16:gp4
  d16

  ; back-transfer
  (p1 ph1)
  4u
  p16:gp5
  d16
  DELTA1 pl0:f2
  (center (p2 ph1) (p8:sp13 ph1):f2 )
  4u
  p16:gp6
  DELTA1 pl12:f2
  DELTA3 BLKGRAD
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 
     F1EA(calgrad(EA), caldel(d0, +in0) & calph(ph12, +180) & calph(ph31, +180))



exit
   

ph1=0 
ph2=1
ph11=0 2
ph12=0 0 2 2
ph13=0 0 0 0 2 2 2 2
ph6=0
;ph30=0
ph31=0 2 0 2 2 0 2 0
  

;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse
;p19: second homospoil/gradient pulse
;p22: f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                  [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(XH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 1 * n
;DS: >= 16
;td1: number of experiments
;FnMODE: echo-antiecho
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:	gp 1 : gp 2
;			  80 : 20.1    for C-13
;			  80 :  8.1    for N-15

;for z-only gradients:
;gpz1: 80%
;gpz2: 20.1% for C-13, 8.1% for N-15
;gpz3: 7%

;use gradient files:   
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: hsqcetgp,v 1.4 2007/04/11 13:34:30 ber Exp $
