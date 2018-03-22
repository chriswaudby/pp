;Delays corrected to give zero first-order phase correction
;With gradients during t1 to keep water along z
;
;hsqcetfpf3gpsi2
;avance-version (07/04/04)
;HSQC
;2D H-1/X correlation via double inept transfer
;   using sensitivity improvement
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;using f3 - channel
;using flip-back pulse
;with gradients in back-inept
;
;A.G. Palmer III, J. Cavanagh, P.E. Wright & M. Rance, J. Magn.
;   Reson. 93, 151-170 (1991)
;L.E. Kay, P. Keifer & T. Saarinen, J. Am. Chem. Soc. 114,
;   10663-5 (1992)
;J. Schleucher, M. Schwendinger, M. Sattler, P. Schmidt, O. Schedletzky,
;   S.J. Glaser, O.W. Sorensen & C. Griesinger, J. Biomol. NMR 4,
;   301-306 (1994)
;S. Grzesiek & A. Bax, J. Am. Chem. Soc. 115, 12593-12594 (1993)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d26=1s/(cnst4*4)"


"d0=6u"

"in0=inf1/4"


"DELTA1=p16+d16+8u+de-0.6366*p1"
"DELTA2=d24-p19-d16"
"DELTA3=d26-p16-d16"
"acqt0=de"

#   ifdef LABEL_CN
"DELTA=p16+d16+larger(p2,p8)+d0*4+20u"
#   else
"DELTA=p16+d16+p2+d0*4+20u"
#   endif /*LABEL_CN*/


1 ze
  d11 pl16:f3
2 d1 do:f3
3 (p1 ph1)
  d26 pl3:f3
  (center (p2 ph2) (p22 ph6):f3 )
  d26 UNBLKGRAD
  (p1 ph2) 
  4u pl0:f1
  (p11:sp1 ph1:r):f1		; flipdown, -y -> -z
  4u
  p16:gp1
  d16 pl1:f1
  (p21 ph3):f3
  
  2u
  d0 gron0
  d0 gron0*-1
  8u groff 

#   ifdef LABEL_CN
  (center (p2 ph7) (p8:sp13 ph1):f2 )
#   else
  (p2 ph7)
#   endif /*LABEL_CN*/

  2u
  d0 gron0
  d0 gron0*-1
  8u groff

  p16:gp2*EA
  d16
  (p22 ph4):f3
  DELTA
  (center (p1 ph1) (p21 ph4):f3 )
  p19:gp4
  d16
  DELTA2
  (center (p2 ph1) (p22 ph1):f3 )
  DELTA2
  p19:gp4
  d16
  (center (p1 ph2) (p21 ph5):f3 )
  p16:gp5
  d16
  DELTA3
  (center (p2 ph1) (p22 ph1):f3 )
  DELTA3
  p16:gp5
  d16
  (p1 ph1)
  DELTA1
  (p2 ph1)
  4u
  p16:gp3
  d16 pl16:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3
  d1 do:f3 mc #0 to 2 
     F1EA(igrad EA & ip5*2, id0 & ip3*2 & ip6*2 & ip31*2)
exit
   

ph1=0 
ph2=1
ph3=0 2
ph4=0 0 2 2
ph5=1 1 3 3
ph6=0
ph7=0 0 2 2
ph31=2 0 0 2
  

;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl16: f3 channel - power level for CPD/BB decoupling
;sp1: f1 channel - shaped pulse  90 degree [flipdown]
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse                         [1 msec]
;p19: gradient pulse 2                                 [500 usec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                           [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery
;d24: 1/(4J)YH for YH
;     1/(8J)YH for all multiplicities
;d26: 1/(4J(YH))
;cnst4: = J(YH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 1 * n
;DS: >= 16
;td1: number of experiments
;FnMODE: echo-antiecho
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:  gp 1 : gp 2 : gp 3 : gp 4 : gp 5
;                       50 :   80 : 20.1 :    5 :   -2       for C-13
;                       50 :   80 :  8.1 :    5 :   -2       for N-15

;for z-only gradients:
;gpz0: 1-2%
;gpz1: 50%
;gpz2: 80%
;gpz3: 8.1% for N-15, 20.1% for C-13
;gpz4: 5%
;gpz5: -2%

;use gradient files:   
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam5: SINE.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: hsqcetfpf3gpsi2,v 1.3 2007/04/11 13:34:30 ber Exp $
