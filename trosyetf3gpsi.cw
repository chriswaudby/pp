;With distinct power levels for flipdown/flipback pulses.
;
;trosyetf3gpsi.2
;avance-version (08/01/15)
;2D H-1/X correlation via TROSY
;   using sensitivity improvement
;phase sensitive using Echo/Antiecho gradient selection
;using f3 - channel
;(use parameterset )
;
;T. Schulte-Herbrueggen & O.W. Sorensen, J. Magn. Reson. 144,
;   123 - 128 (2000)
;(M. Czisch & R. Boelens, J. Magn. Reson. 134, 158-160 (1998) )
;(K. Pervushin, G. Wider & K. Wuethrich, J. Biomol. NMR 12,
;   345-348 (1998) )
;(A. Meissner, T. Schulte-Herbrueggen, J. Briand & O.W. Sorensen, Mol. Phys. 96,
;   1137-1142 (1998) )
;(J. Weigelt, J. Am. Chem. Soc. 120, 10778-10779 (1998) )
;(M. Rance, J.P. Loria & A.G. Palmer III, J. Magn. Reson. 136, 91-101 (1999) )
;(G. Zhu, X.M. Kong & K.H. Sze, J. Biomol. NMR 13, 77-81 (1999) )
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


define list<gradient> EA2 = { 0.8750 1.0000}
define list<gradient> EA4 = { 1.0000 0.6667}
define list<gradient> EA6 = { 0.6595 1.0000}


"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"
"d13=4u"
"d26=1s/(cnst4*4)"


"in0=inf1/2"

"d0=20u"


"DELTA1=d26-p16-d16"
"DELTA2=d25-p16-d16"
"DELTA3=d26-p11-p16-d16-8u"

#   ifdef LABEL_CN
"DELTA=d0*2+p8-p21*4/3.1416+8u"
#   else
;"DELTA=d0*2-p21*4/3.1416+6u"
"DELTA=d0*2+6u"
#   endif /*LABEL_CN*/
"acqt0=0"

1 ze 
2 d11
3 d1 pl1:f1
  50u UNBLKGRAD

  (p1 ph3)
  p16:gp1
  d16
  DELTA1
  (center (p2 ph2) (p22 ph1):f3 )
  DELTA1
  p16:gp1
  d16
  (p1 ph2) 

  p16:gp2*EA2
  d16

  (p21 ph5):f3
  DELTA
  (p22 ph1):f3

#   ifdef LABEL_CN
  d0 gron0
  2u groff
  (p8:sp13 ph1):f2
  d0 gron0*-1
  2u groff
#   else
  d0 gron0
  d0 gron0*-1
  2u groff
#   endif /*LABEL_CN*/

  4u
  p16:gp2*-1*EA2
  d16

  (p1 ph6)
  p16:gp3
  d16
  DELTA2
  (center (p2 ph2) (p22 ph2):f3 )
  DELTA2
  p16:gp3
  d16
  (p1 ph1)

  p16:gp4*EA4
  d16

  (p21 ph7):f3
  p16:gp5
  d16
  DELTA3 pl0:f1
  (p11:sp1 ph4:r):f1	; flipdown(-y), z -> -x
  4u
  4u pl1:f1
  (center (p2 ph2) (p22 ph8):f3 )
  4u pl0:f1
  (p11:sp11 ph14:r):f1	; flipback(-y), x -> z
  4u
  DELTA3
  p16:gp5
  d16 pl1:f1
  (p21 ph9:r):f3

  p16:gp6*EA6
  d16
  4u BLKGRAD

  go=2 ph31
  d11 mc #0 to 2 
     F1EA(igrad EA2 & igrad EA4 & igrad EA6 & ip6*2 & ip7*2, id0 & ip5*2 & ip31*2)
exit 
  

ph1=0
ph2=1
ph3=2
ph4=3
ph14=3
ph5=0 2
ph6=3
ph7=0 0 2 2
ph8=1 1 3 3
ph9=3 3 1 1
ph31=0 2 2 0


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;sp1: f1 channel - shaped pulse  90 degree  (flipdown)
;sp11: f1 channel - shaped pulse  90 degree  (flipback)
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p11: f1 channel -  90 degree shaped pulse             [1 msec]
;p16: homospoil/gradient pulse                         [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                           [6 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery
;d25 : 1/(4J')NH, 
;         compensation delay for suppression of other Trosy peaks
;d26 : 1/(4J)NH
;cnst4: = J(NH)
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(N)) = DW(N)
;nd0: 2
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: echo-antiecho


;use gradient ratio: gp 0 : gp 1 : gp 2 : gp 3 : gp 4 : gp 5 : gp 6
;                       3 :    3 :   80 :    5 :   30 :    7 :30.13

;for z-only gradients:
;gpz0: 3%
;gpz1: 3%
;gpz2: 80%
;gpz3: 5%
;gpz4: 30%
;gpz5: 7%
;gpz6: 30.13%

;use gradient files:
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam5: SINE.100
;gpnam6: SINE.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;Processing

;PHC0(F1): 22.5



;$Id: trosyetf3gpsi.2,v 1.1.2.1 2008/01/15 18:13:12 ber Exp $
