;sfhmqcf3gpph
;avance-version (16/01/25)
;SOFAST HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive
;with decoupling during acquisition
;use baseopt
;
;P.Schanda and B. Brutscher, J. Am. Chem. Soc. 127, 8014 (2005)
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


"d11=30m"
"d12=20u"
"d21=1s/(cnst4*2)"


#   ifdef CALC_SP
"p39=(bwfac23/(cnst55*cnst49*bf1))*1000000"
"spw23=plw1/((p39*90.0)/(p1*totrot23))*((p39*90.0)/(p1*totrot23))*(integfac23*integfac23)"
"spoal23=1"

"p40=(bwfac24/(cnst55*cnst50*bf1))*1000000"
"spw24=plw1/((p40*90.0)/(p1*totrot24))*((p40*90.0)/(p1*totrot24))*(integfac24*integfac24)"
"spoal24=0.5"
#   endif /*CALC_SP*/


"in0=inf1"

"d0=in0/2-p21*4/3.1415"


"DELTA1=d21-p16-d16-p39*cnst39"
"DELTA2=p39*cnst39-4u"
"acqt0=0"

"spoff23=bf1*(cnst54/1000000)-o1"
"spoff24=bf1*(cnst54/1000000)-o1"


1 ze 
  d11 pl26:f3
2 d1 do:f3
3 d12 pl3:f3
  50u UNBLKGRAD

  p16:gp2
  d16

  (p39:sp23 ph1):f1
  p16:gp1
  d16

#   ifdef LABEL_CN
  (center (p40:sp24 ph2):f1 (p8:sp13 ph1):f2 (DELTA1 p21 ph3 d0 p21 ph4 DELTA1):f3 )
#   else
  (center (p40:sp24 ph2):f1 (DELTA1 p21 ph3 d0 p21 ph4 DELTA1):f3 )
#   endif /*LABEL_CN*/


  DELTA2
  p16:gp1
  d16 pl26:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3 
  d1 do:f3 mc #0 to 2 
     F1PH(calph(ph3, +90), caldel(d0, +in0))
exit 
  

ph1=0 
ph2=0 
ph3=0 2
ph4=0 0 2 2 
ph31=0 2 2 0


;pl3 : f3 channel - power level for pulse (default)
;pl26: f3 channel - power level for CPD/BB decoupling (low power)
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;sp23: f1 channel - shaped pulse 120 degree 
;                   (Pc9_4_120.1000)
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse                       [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (2.4ms at 600.13 MHz)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000               (0.8ms at 600.13 MHz)
;d0 : incremented delay (2D) = in0/2-p21*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)NH
;cnst4: = J(NH)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;cnst49: scaling factor for p39 to compensate for transition region
;           Pc9_4_120.1000: 1.074
;cnst50: scaling factor for p40 to compensate for transition region
;           Rsnob.1000: 1.000
;cnst54: H(N) chemical shift (offset, in ppm)
;cnst55: H(N) bandwidth (in ppm)
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/ SW(N) = 2 * DW(N)
;nd0: 1
;ns: 2 * n
;ds: 16
;aq: <= 50 msec
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEC
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence
;          use pulse of >= 350 usec


;use gradient ratio:	gp 1 : gp 2
;			  11 :    7


;for z-only gradients:
;gpz1: 11%
;gpz2:  7%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


                                          ;preprocessor-flags-start
;CALC_SP: for calculation of all bandselective Proton pulses based on cnst54 and cnst55
;             option -DCALC_SP (eda: ZGOPTNS)
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1



;$Id: sfhmqcf3gpph,v 1.12.2.4 2016/04/06 09:21:25 ber Exp $
