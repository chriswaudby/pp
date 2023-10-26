;project_hsqcetgpsp.2
;13C gHSQC incorporating a PROJECT CPMG sequence for 1H T2 relaxation measurement
;
;avance-version (12/01/11)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;using trim pulses in inept transfer
;using shaped pulses for inversion and refocussing on f2 - channel
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"d4=1s/(cnst2*4)"
"d11=30m"

#   ifdef LABEL_CN
"p22=p21*2"
#   else
#   endif /*LABEL_CN*/


"d0=3u"

"in0=inf2/2"

define list<loopcounter> ncyc=<$VCLIST>

"DELTA1=d4-p16-larger(p2,p14)/2-de-8u"
"DELTA2=d4-larger(p2,p14)/2"
"DELTA3=d4-larger(p2,p14)/2-p1*2/PI"

"d2=0.5ms"
"d20=d2-0.5*p4"
"d21=d2-0.5*p4-p1"
"d22=d2-0.5*p4-0.5*p1"


#   ifdef LABEL_CN
"DELTA=p16+d16+larger(p2,p22)+d0*2"
#   else
"DELTA=p16+d16+p2+d0*2"
#   endif /*LABEL_CN*/

"l1=0"
"l2=0"

"acqt0=0"
baseopt_echo
aqseq 312

1 ze
  d11 pl12:f2
2 d11 do:f2

  ; off-res presat
  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  d1 cw:f1 ph1
  4u do:f1
  30u fq=0:f1

  "l2 = (trunc(ncyc[l1] + 0.3))"  
  4u pl1:f1
  4u pl2:f2

  (p1 ph1)

  ; PROJECT
3 d20
  (p4 ph1):f2
  d21
  p2 ph2
  d21
  (p4 ph1):f2
  d22
  p1 ph7
  d22
  (p4 ph1):f2
  d21
  p2 ph2
  d21
  (p4 ph1):f2
  d20
  lo to 3 times l2
 
  ; continue HSQC
  DELTA2 pl0:f2
  4u
  (center (p2 ph1) (p14:sp3 ph6):f2 )
  4u
  DELTA2 pl2:f2 UNBLKGRAD
  p28 ph1
  4u
  (p1 ph2) (p3 ph3):f2
  d0 

#   ifdef LABEL_CN
  (center (p2 ph5) (p22 ph1):f3 )
#   else
  (p2 ph5)
#   endif /*LABEL_CN*/

  d0
  p16:gp1*EA
  d16 pl0:f2
  4u
  (p24:sp7 ph4):f2
  4u
  DELTA pl2:f2
  (ralign (p1 ph1) (p3 ph4):f2 )
  DELTA3 pl0:f2
  (center (p2 ph1) (p14:sp3 ph1):f2 )
  4u
  p16:gp2
  DELTA1 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2
  d11 do:f2 mc #0 to 2 
     F1QF(calclc(l1,1))
     F2EA(calgrad(EA), caldel(d0, +in0) & calph(ph3, +180) & calph(ph6, +180) & calph(ph31, +180))
exit
   

ph1=0 
ph2=1
ph3=0 2
ph4=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph5=0 0 2 2
ph6=0
ph7=1 1 1 1 3 3 3 3
ph31=0 2 0 2 0 2 0 2 2 0 2 0 2 0 2 0


;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;sp3: f2 channel - shaped pulse 180 degree for inversion
;sp7: f2 channel - shaped pulse 180 degree for refocussing
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p14: f2 channel - 180 degree shaped pulse for inversion
;p16: homospoil/gradient pulse
;p22: f3 channel - 180 degree high power pulse
;p24: f2 channel - 180 degree shaped pulse for refocussing
;p28: f1 channel - trim pulse
;d0 : incremented delay (2D)                         [3 usec]
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



;$Id: hsqcetgpsp.2,v 1.9 2012/01/31 17:49:26 ber Exp $

