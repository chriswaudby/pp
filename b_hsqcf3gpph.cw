;BEST-15N HSQC, using PC-9 shapes for all 90deg pulses on Hn
;Chris W, Aug 2018
;
;b_hsqcf3gpph.cw
;avance-version (08/01/24)
;(best)-HSQC
;2D H-1/X correlation via double inept transfer
;   using sensitivity improvement
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;using f3 - channel
;with gradients in back-inept
;
;(A.G. Palmer III, J. Cavanagh, P.E. Wright & M. Rance, J. Magn.
;   Reson. 93, 151-170 (1991) )
;(L.E. Kay, P. Keifer & T. Saarinen, J. Am. Chem. Soc. 114,
;   10663-5 (1992) )
;(J. Schleucher, M. Schwendinger, M. Sattler, P. Schmidt,
;   O. Schedletzky, S.J. Glaser, O.W. Sorensen & C. Griesinger, 
;   J. Biomol. NMR 4, 301-306 (1994) )
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


"p22=p21*2"
"d11=30m"
"d12=20u"
"d26=1s/(cnst4*4)"

"d0=inf1/2-2*0.6366*p21"
"in0=inf1"

"DELTA1=d26-p19-d16-p41*cnst41-p42*cnst42"
"DELTA2=DELTA1+p41*cnst41-4u-de"

"spoff25=bf1*(cnst19/1000000)-o1"
"spoff26=bf1*(cnst19/1000000)-o1"
"spoff27=bf1*(cnst19/1000000)-o1"
"d42=p42"

"acqt0=de"


1 ze
  d11 pl26:f3 pl0:f1
2 d11 do:f3
  d1
  50u UNBLKGRAD
  4u pl3:f3

  ; first INEPT
  (p41:sp25 ph1)
  DELTA1
  p19:gp1
  d16
  (center (p42:sp26 ph1) (p22 ph1):f3 )
  DELTA1
  p19:gp1
  d16
  (p41:sp27 ph2):f1

  ; zz filter
  p16:gp2
  d16

  ; t1 evolution
#   ifdef LABEL_CN
  (center
    (p42:sp26 ph5):f1
    (p8:sp13 ph1):f2
    (d42*0.5 p21 ph3 d0 p21 ph4 d42*0.5):f3
  )
#   else
  (center
    (p42:sp26 ph5):f1
    (d42*0.5 p21 ph3 d0 p21 ph4 d42*0.5):f3
  )
#   endif /*LABEL_CN*/

  ; zz filter
  p16:gp3
  d16

  ; retro-INEPT
  (p41:sp25 ph1)
  DELTA1
  p19:gp1
  d16
  (center (p42:sp26 ph1) (p22 ph1):f3 )
  DELTA2 pl26:f3
  p19:gp1
  d16
  4u BLKGRAD 

  ; acquisition
  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
     F1PH(ip3, id0)
exit
   

ph1=0 
ph2=1
ph3=0 2
ph4=0 0 2 2
ph5=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph31=0 2 2 0
  

;pl3 : f3 channel - power level for pulse (default)
;pl26: f3 channel - power level for CPD/BB low power decoupling
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;sp25: f1 channel - shaped pulse  90 degree (Pc9_4_90.1000)
;sp26: f1 channel - shaped pulse 180 degree (Reburp.1000)
;sp27: f1 channel - shaped pulse  90 degree (Pc9_4_90.1000)
;                   for time reversed pulse
;sp30: f1 channel - shaped pulse 180 degree (Bip720,50,20.1)
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse                         [1 msec]
;p19: gradient pulse 2                                 [500 usec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p29: gradient pulse 3                                 [250 usec]
;p41: f1 channel -  90 degree shaped pulse for excitation
;                      Pc9_4_90.1000             (2230 usec at 600 MHz)
;p42: f1 channel - 180 degree shaped pulse for refocusing
;                      Reburp.1000               (1680 usec at 600 MHz)
;p43: f1 channel -  90 degree shaped pulse for excitation
;                      Bip720,50,20.1            (200us at 600 MHz)
;d0 : incremented delay (2D)                           [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery
;d26: 1/(4J(NH))
;d63: set to zero unless DELTA becomes negative
;cnst4: = J(NH)
;cnst19: H(N) chemical shift (offset, in ppm)
;cnst41: compensation of coupling evolution during p41
;           Pc9_4_90.1000: 0.514
;cnst42: compensation of coupling evolution during p42
;           Reburp.1000: 0.475
;NS: 2 * n
;DS: >= 16
;aq: <= 50 msec (or <100ms with d1>100ms)
;td1: number of experiments
;FnMODE: echo-antiecho
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence

;for z-only gradients:
;gpz1: 8%
;gpz2: 13% for N-15
;gpz3: -7%

;use gradient files:   
;gpnam1: SMSQ10.50
;gpnam2: SMSQ10.50
;gpnam3: SMSQ10.50

                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end

