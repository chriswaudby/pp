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

"DELTA1=d4-p16-d16-4u-larger(p2,p4)*0.5"
"DELTA2=2*d0+p2"
"DELTA3=d16+0.6366*p1"
"acqt0=0"

1 ze
  d11 pl1:f1 pl12:f2
2 d11 do:f2

  d1
  50u UNBLKGRAD
  4u pl1:f1 pl2:f2
 
   ; purge Cz
  (p3 ph1):f2
  4u
  p16:gp0
  d16

   ; begin HSQC
  (p1 ph1)
  4u
  p16:gp1
  d16
  DELTA1
  (center (p2 ph2) (p4 ph1):f2 )
  DELTA1 
  p16:gp1
  d16 
  4u
  (p1 ph2)

#ifdef CALIB_C
  ; zz filter
  4u
  p16:gp2
  d16

  4u pl20:f2
  (p20 ph1):f2
  4u pl2:f2
#endif

#ifdef ZZ1
  ; zz filter
  4u
  p16:gp2
  d16
#endif

  ; t1 evolution
  (p3 ph11):f2
  d0
  (p2 ph13)
  d0
  p16:gp3*EA
  d16
  (p4 ph14):f2
  p16:gp3*EA*-1
  d16
  DELTA2
  (p3 ph12):f2

#ifdef ZZ2
  ; zz filter
  4u
  p16:gp4
  d16
#endif

  ; back-transfer
  (p1 ph1)
  4u
  p16:gp5
  d16
  DELTA1 
  (center (p2 ph1) (p4 ph1):f2 )
  4u
  p16:gp6
  DELTA1 pl12:f2
  DELTA3 BLKGRAD
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 
     F1EA(igrad EA, id0)



exit
   

ph1=0 
ph2=1
ph11=0 2
ph12=0 0 2 2
ph13=0 0 0 0 2 2 2 2
ph14=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph6=0
ph31=0 2 2 0
  

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

;for z-only gradients:
;gpz1: 7%
;gpz2: 19%
;gpz3: 80%
;gpz4: -23%
;gpz5: 30.1%
;gpz6: 10%

;use gradient files:   
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam5: SINE.100
;gpnam6: SINE.100

