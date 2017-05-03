;calib13C.cw
;based on hsqcetgp
;avance-version (12/01/11)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;using trim pulses in inept transfer
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p4=p3*2"
"d4=1s/(cnst2*4)"
"d11=30m"

"d0=3u"



"DELTA1=d4-p16-de+p1*2/PI-8u"


"DELTA=p16+d16+p2+d0*2"



1 ze
  d11 pl1:f1 pl12:f2
2 d1 do:f2 
3 (p1 ph1)
  d4 pl2:f2
  (center (p2 ph1) (p4 ph6):f2 )
  d4 UNBLKGRAD
  (p1 ph2)
p16:gp3
d16
# ifdef CAL_C
4u pl20:f2
(p20 ph1):f2
p16:gp3*1.33
d16 pl2:f2
# endif /*CAL_C*/

(p3 ph3):f2
  3u
  (p2 ph5)
  3u
  p16:gp1*EA
  d16
  (p4 ph4):f2
  DELTA
  (ralign (p1 ph1) (p3 ph4):f2 )
  d4
  (center (p2 ph1) (p4 ph1):f2 )
  4u
  p16:gp2
  DELTA1 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2
  d1 do:f2 mc #0 to 2 
exit
   

ph1=0 
ph2=1
ph3=0 2
ph4=0 0 0 0 2 2 2 2
ph5=0 0 2 2
ph6=0
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
;p22: f3 channel - 180 degree high power pulse
;p28: f1 channel - trim pulse
;d0 : incremented delay (2D)                  [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(XH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;ns: 1 * n
;ds: >= 16
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

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: hsqcetgp,v 1.5.4.1.4.1 2012/01/31 17:56:32 ber Exp $
