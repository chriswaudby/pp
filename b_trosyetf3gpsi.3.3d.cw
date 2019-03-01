;b_trosyetf3gpsi.3
;pseudo-3D for kinetics
;
;avance-version (15/03/12)
;best-TROSY
;2D H-1/X correlation via TROSY
;   using sensitivity improvement
;phase sensitive using Echo/Antiecho 
;using f3 - channel
;using shaped pulses for inversion and refocussing on f3
;uncompensated version d25=d26
;with additional 180degree pulse on N-15
;(use parameterset B_TROSYETF3GPSI)
;
;Z. Solyom, M. Schwarten, L. Geist, R. Konrat D. Willbold &
;   Bernhard Brutscher, J. Biomol. NMR 55, 311-321 (2013)
;A. Favier & B. Brutscher, J. Biomol. NMR 49, 9-15 (2011)
;(E. Lescop, P. Schanda & B. Brutscher,
;   J. Magn. Reson.  187 163-169 (2007))
;(T. Schulte-Herbrueggen & O.W. Sorensen, J. Magn. Reson. 144, 
;   123 - 128 (2000))
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


define list<gradient> EA3 = { 1.0000 0.8750 }
define list<gradient> EA5 = { 0.6667 1.0000 }
define list<gradient> EA7 = { 1.0000 0.6595 }


"d11=30m"

"d25=2.7m"
"d26=2.7m"

"p19=500u"
"p29=250u"


#   ifdef CALC_SP
"p41=(bwfac25/(cnst55*cnst51*bf1))*1000000"
"spw25=plw1/((p41*90.0)/(p1*totrot25))*((p41*90.0)/(p1*totrot25))*(integfac25*integfac25)"
"spw27=plw1/((p41*90.0)/(p1*totrot27))*((p41*90.0)/(p1*totrot27))*(integfac27*integfac27)"
"spoal25=1"
"spoal27=0"

"p42=(bwfac26/(cnst55*cnst52*bf1))*1000000"
"spw26=plw1/((p42*90.0)/(p1*totrot26))*((p42*90.0)/(p1*totrot26))*(integfac26*integfac26)"
"spoal26=0.5"

"p43=(bwfac28/(cnst55*cnst53*bf1))*1000000"
"spw28=plw1/((p43*90.0)/(p1*totrot28))*((p43*90.0)/(p1*totrot28))*(integfac28*integfac28)"
"spw29=plw1/((p43*90.0)/(p1*totrot29))*((p43*90.0)/(p1*totrot29))*(integfac29*integfac29)"
"spoal28=1"
"spoal29=0"
#   endif /*CALC_SP*/


"d0=3u"

"in0=inf1/2"


"DELTA1=d26-p29-d16-larger(p56,p42)/2-p41*cnst41"
"DELTA6=d25-p29-d16-larger(p56,p42)/2-p43*cnst43"
"DELTA7=d26-p16-d16-larger(p57,p42)/2"
"DELTA8=de+4u"

#   ifdef LABEL_CN
"DELTA=d0*2+p8+p21*4/PI+de+4u"
#   else
"DELTA=d0*2+p21*4/PI+de+4u"
#   endif /*LABEL_CN*/


"spoff13=bf2*(cnst26/1000000)-o2"

"spoff25=bf1*(cnst54/1000000)-o1"
"spoff26=bf1*(cnst54/1000000)-o1"
"spoff28=bf1*(cnst54/1000000)-o1"
"spoff27=bf1*(cnst54/1000000)-o1"
"spoff29=bf1*(cnst54/1000000)-o1"

"spoff39=0"
"spoff40=0"


"acqt0=0"
baseopt_echo

aqseq 321

1 d11 ze
2 3m
  
  (p56:sp39 ph1):f3
  d1
  50u UNBLKGRAD

  (p41:sp25 ph1)
  p29:gp1
  d16
  DELTA1
  (center (p42:sp26 ph1) (p56:sp39 ph1):f3 )
  DELTA1
  p29:gp1
  d16
  (p41:sp27 ph2):f1

  p29:gp2
  d16 pl3:f3

  (p21 ph5):f3
  d0

#   ifdef LABEL_CN
  (p8:sp13 ph1):f2
#   else
#   endif /*LABEL_CN*/

  d0
  (p56:sp39 ph1):f3
  DELTA

  p19:gp3*EA3
  d16

  (p43:sp29 ph6)
  p29:gp4
  d16
  DELTA6
  (center (p42:sp26 ph1) (p56:sp39 ph1):f3 )
  DELTA6
  p29:gp4
  d16
  (p43:sp28 ph1)

  p19:gp5*EA5
  d16 pl3:f3
  DELTA8

  (p21 ph2):f3
  p16:gp6
  d16
  DELTA7
  (center (p42:sp26 ph1) (p57:sp40 ph1):f3 )
  DELTA7
  p16:gp6
  d16 pl3:f3
  (p21 ph7:r):f3

  p19:gp7*EA7
  d16
  4u BLKGRAD

  go=2 ph31
  3m mc #0 to 2 
     F2EA(calgrad(EA3) & calgrad(EA5) & calgrad(EA7) & calph(ph6, +180) & calph(ph7, +180), caldel(d0, +in0) & calph(ph5, +180) & calph(ph31, +180))
     F1QF()
exit


ph1=0
ph2=1 
ph3=2
ph4=3
ph5=0 2
ph6=3
ph7=2
ph31=0 2


;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;sp13: f2 channel - shaped pulse 180 degree (Ca and C=O, adiabatic)
;sp25: f1 channel - shaped pulse  90 degree (Pc9_4_90.1000)
;sp26: f1 channel - shaped pulse 180 degree (Reburp.1000)
;sp27: f1 channel - shaped pulse  90 degree (Pc9_4_90.1000)
;                   for time reversed pulse
;sp28: f1 channel - shaped pulse  90 degree (Eburp2.1000)
;sp29: f1 channel - shaped pulse  90 degree (Eburp2tr.1000)
;                   for time reversed pulse
;sp39: f3 channel - shaped pulse 180 degree (Bip720,50,20.1)
;sp40: f3 channel - shaped pulse 180 degree (Reburp.1000)
;sp26: f1 channel - shaped pulse 180 degree (Reburp.1000)
;sp28: f1 channel - shaped pulse  90 degree (Eburp2.1000)
;sp29: f1 channel - shaped pulse  90 degree (Eburp2tr.1000)
;                   for time reversed pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse                         [1 msec]
;p19: gradient pulse 2                                 [500 usec]
;p21: f3 channel -  90 degree high power pulse
;p29: gradient pulse 3                                 [250 usec]
;p41: f1 channel -  90 degree shaped pulse for excitation
;                      Pc9_4_90.1000             (2.2ms at 600.13 MHz)
;p42: f1 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000               (1.4ms at 600.13 MHz)
;p43: f1 channel -  90 degree shaped pulse for excitation
;                      Eburp2.1000/Eburp2tr.1000 (1.7ms at 600.13 MHz)
;p56: f3 channel - 180 degree shaped pulse for inversion
;                      Bip720,50,20.1            (500us at 600.13 MHz)
;p57: f3 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000               (1.6ms at 600.13 MHz)
;p42: f1 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000               (2.0ms at 600.13 MHz)
;p43: f1 channel -  90 degree shaped pulse for excitation
;                      Eburp2.1000/Eburp2tr.1000 (1.92ms at 600.13 MHz)
;p56: f3 channel - 180 degree shaped pulse for inversion
;p57: f3 channel - 180 degree shaped pulse for refocussing
;d0 : incremented delay (F1)                           [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery
;d25: 1/(4J'(NH)                                       [2.7 msec]
;d26: 1/(4J(NH)                                        [2.7 msec]
;cnst26: Call chemical shift (offset, in ppm)          [101 ppm]
;cnst41: compensation of chemical shift evolution during p41
;           Pc9_4_90.1000: 0.529
;cnst43: compensation of chemical shift evolution during p43
;           Eburp2.1000: 0.69
;cnst51: scaling factor for p41 to compensate for transition region
;           Pc9_4_90.1000: 1.172
;cnst52: scaling factor for p42 to compensate for transition region
;           Reburp.1000: 1.426
;cnst53: scaling factor for p43 to compensate for transition region
;           Eburp2.1000: 1.000
;cnst54: H(N) chemical shift (offset, in ppm)
;cnst55: H(N) bandwidth (in ppm)
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(N)) = DW(N)
;nd0: 2
;ns: 2 * n
;ds: >= 16
;td1: number of experiments
;FnMODE: echo-antiecho


;for z-only gradients:
;gpz1: 2%
;gpz2: 21%
;gpz3: -80%
;gpz4: 5%
;gpz5: 30%
;gpz6: 45%
;gpz7: 30.13%

;use gradient files:   
;gpnam1: SMSQ10.32
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.32
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam7: SMSQ10.100



                                          ;preprocessor-flags-start
;CALC_SP: for calculation of all bandselective Proton pulses based on cnst54 and cnst55
;             option -DCALC_SP (eda: ZGOPTNS)
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end
										  

										  
;Processing

;PHC0(F1): 45.0


										  
;$Id: b_trosyetf3gpsi.3,v 1.1.2.2 2015/03/12 17:07:08 ber Exp $
