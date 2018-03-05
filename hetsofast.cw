;BB_2d_hetsofast_hmqc
;correct version - run this one!
;well-calibrated Reburp (1942.8 us @ 700 MHz, 8.2ppm)
;and aliphatic inversion (2570 us @ 700 MHz, 1.7ppm)
;allowing short recycle delays for aliphatic and reference spectra
;
;Chris Waudby, Anais Cassaignau Aug 2015
;
;$CLASS=IBS
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=

prosol relations=<IBS>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define delay dNH

"d11=30m"
"d12=20u"
"dNH=1s/(cnst4*2)"


"in0=inf1"

"d0=in0/2-p21*4/3.1415"
/*******************************************************************/
/*   calculation of shaped 1H pulse parameters                     */
/*******************************************************************/

"p41=7.2/(cnst2*bf1/1000000)" /*  PC9  pulse length  */
"spw25=plw1*(pow((cnst3/90.0)*(p1/p41)/0.125,2))" /* PC9  power level  */
"spoff25=bf1*(cnst1/1000000)-o1"  /*  PC9  offset */
;"spoal25=0.5"

"p42=1.08*4.875/(cnst2*bf1/1000000)" /* REBURP pulse length  */
"spw26=plw1*(pow((p1*2/p42)/0.0798,2))"   /* REBURP power level  */
"spoff26=bf1*(cnst1/1000000)-o1" /* REBURP offset */
;"spoal26=0.5"

"p45=1.6*4.5/(4.0*bf1/1000000)" /* ISNOB5 pulse length  */
"spw31=plw1*(pow((p1*2.0/p45)/0.0914,2))"   /* ISNOB5 power level  */
"spoff31=bf1*(cnst5/1000000)-o1" /* ISNOB5 offset */
"spoal31=0.5"
 
 
"p46=150000*(600/bf1)" /* ISNOB5 pulse length  */
"spw32=plw1*(pow((p1*2.0/p46)/0.0914,2))"   /* ISNOB5 power level  */
"spoff32=0" /* ISNOB5 offset */
"spoal32=0.5"

/*******************************************************************/
/*   13C Adiabatic pulse                                           */
/*******************************************************************/
#ifdef LABEL_CN
"p8=500u"       /*      SPNAM 13","Crp80,0.5,20.1   */
"if ( bf2 < 165 ) {spw13=plw2*(pow((p3/25.5832),2)); } else{spw13=plw2*(pow((p3/22.1557),2));} "
#endif /*LABEL_CN*/

/*******************************************************************/
/*   calculation of 15N decoupling                                 */
/*******************************************************************/
;"pcpd3=300u"
;"plw26 =plw3*(pow((p21/pcpd3),2))" /* CPD  power level  */
"plw26 =plw3*(pow((p21/p62),2))" /* CPD  power level  */

"DELTA1=dNH-p16-d16-p41*cnst39"
"DELTA2=p41*cnst39-de-4u"
"DELTA3=d1-p45"
"if ( d1 > p46 ) { d23=d1-p46; } else { d23 = p46; }"
"DELTA4=d23"




1 ze 
  d11 pl26:f3
2 d11 do:f3
  d12 pl3:f3
/*******************************************/
/*  Options for HETSOFAST                  */
/*******************************************/
 if "l1==1" 
 {
 d1 
 }
 if "l1==2" 
 {
  DELTA3
  (p45:sp31 ph1):f1
 }
 if "l1==3" 
 {
  (p46:sp32 ph1):f1
  DELTA4
 }
/*******************************************/
  50u UNBLKGRAD

  p16:gp2
  d16

  (p41:sp25 ph1):f1
  p16:gp1
  d16

#   ifdef LABEL_CN
  (center (p42:sp26 ph2):f1 (p8:sp13 ph1):f2 (DELTA1 p21 ph3 d0 p21 ph4 DELTA1):f3 )
#   else
  (center (p42:sp26 ph2):f1 (DELTA1 p21 ph3 d0 p21 ph4 DELTA1):f3 )
#   endif /*LABEL_CN*/


  DELTA2
  p16:gp1
  d16 pl26:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3 
  d11 do:f3 mc #0 to 2 
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



;sp25: f1 channel - shaped pulse 120 degree 
;                   (Pc9_4_120.1000)
;cnst39: compensation of chemical shift evolution during p41/sp25
;           Pc9_4_120.1000: 0.529
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse                       [1 msec]
;p21: 15N -  90 degree high power pulse
;d0 : incremented delay (2D) = in0/2-p21*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;dNH : 1/(2J)NH
;cnst4: = J(NH)
;cnst1: H(N) excitation frequency (in ppm)
;cnst2: H(N) excitation band width (in ppm)
;cnst3:  PC9 flip angle
;cnst5: aliphatic saturation frequency (in ppm)  [1.7 ppm]
;sp26: Reburp.1000
;l1: HETSOFAST_FLG: ref(1) water_sat(3) aliph_sat(2)

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

;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1



