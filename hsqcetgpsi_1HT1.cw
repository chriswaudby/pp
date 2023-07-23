;hsqcetgpsi_1HT1
;avance-version (20/03/03)
;HSQC with 1H inversion recovery
;2D H-1/X correlation via double inept transfer
;   using sensitivity improvement
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;using trim pulses in inept transfer
;
;A.G. Palmer III, J. Cavanagh, P.E. Wright & M. Rance, J. Magn.
;   Reson. 93, 151-170 (1991)
;L.E. Kay, P. Keifer & T. Saarinen, J. Am. Chem. Soc. 114,
;   10663-5 (1992)
;J. Schleucher, M. Schwendinger, M. Sattler, P. Schmidt, O. Schedletzky,
;   S.J. Glaser, O.W. Sorensen & C. Griesinger, J. Biomol. NMR 4,
;   301-306 (1994)
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
"p4=p3*2"
"d4=1s/(cnst2*4)"
"d11=30m"

#   ifdef LABEL_CN
"p22=p21*2"
#   else
#   endif /*LABEL_CN*/


"d0=3u"

"in0=inf2/2"


"DELTA1=p16+d16-p1*0.78+de+8u"

#   ifdef LABEL_CN
"DELTA=p16+d16+50u+larger(p2,p22)+d0*2"
#   else
"DELTA=p16+d16+50u+p2+d0*2"
#   endif /*LABEL_CN*/

define list<delay> t1delay=<$VDLIST>
define delay t1min
"t1min=4*p16+4*d16+4*p3"

"acqt0=0"
baseopt_echo
aqseq 312

1 ze
  t1min
  d11 pl12:f2
2 d1 do:f2 

  "d20=0.25*t1delay-p16-d16-p3"
  "d21=0.25*t1delay-p3"
  4u UNBLKGRAD
  ; 1H inversion (leave water on +z)
  (p10:sp1 ph1):f1
  4u pl1:f1 pl2:f2
  p2 ph1
  p16:gp0
  d16
  ; inversion-recovery period with 13C decoupling to suppress cross-correlated relaxation
  d20 BLKGRAD
  (p4 ph1):f2
  d21
  d21
  (p4 ph1):f2
  d21
  
  (p1 ph1)
  d4
  (center (p2 ph1) (p4 ph6):f2 )
  d4

#   ifdef TRIMP
  p28 ph1
#   endif /* TRIMP */

  (p1 ph2) (p3 ph3):f2
  d0 

#   ifdef LABEL_CN
  (center (p2 ph7) (p22 ph1):f3 )
#   else
  (p2 ph7)
#   endif /*LABEL_CN*/

  d0
  50u UNBLKGRAD
  p16:gp1*EA
  d16
  (p4 ph4):f2
  DELTA
  (center (p1 ph1) (p3 ph4):f2 )
  d24
  (center (p2 ph1) (p4 ph1):f2 )
  d24
  (center (p1 ph2) (p3 ph5):f2 )
  d4
  (center (p2 ph1) (p4 ph1):f2 )
  d4
  (p1 ph1)
  DELTA1
  (p2 ph1)
  4u
  p16:gp2
  d16 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2
  d1 do:f2 mc #0 to 2 
     F1QF(exec(t1delay.inc))
     F2EA(calgrad(EA) & calph(ph5, +180), caldel(d0, +in0) & calph(ph3, +180) & calph(ph6, +180) & calph(ph31, +180))
exit
   

ph1=0 
ph2=1
ph3=0 2
ph4=0 0 2 2
ph5=1 1 3 3
ph6=0
ph7=0 0 2 2
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
;p22: f3 channel - 180 degree high power pulse
;p28: f1 channel - trim pulse
;d0 : incremented delay (2D)                  [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d16: delay for homospoil/gradient recovery
;d24: 1/(8J)XH for all multiplicities
;     1/(4J)XH for XH
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


;for z-only gradients:
;gpz1: 80%
;gpz2: 20.1% for C-13, 8.1% for N-15

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
;TRIMP: to use trimpulse p28@pl1 start experiment with
;          option -DTRIMP (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id:$
