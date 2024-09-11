;3D_NOE_CHSQC_BB
;based on SE 1H-13C HSQC
;with CO decoupling
;with presaturation
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=

prosol relations=<IBS>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p4=p3*2"
"d4=1s/(cnst2*4)"
"d24=1s/(cnst2*8)"
"d11=30m"



/*******************************************************************/
/*   calculation of shaped 13C pulse parameters                    */
/*******************************************************************/
"cnst22 = (sfo2-bf2)*1000000/bf2"    /*  CACB  frequency offset   */
"cnst21 = 173.0"     /*  CO frequency offset   */

"p25=300u"    /* BIP pulse length  */
"spw25=plw2*(pow((p3*8/p14),2))"   /* BIP power level  */
"spoffs25=0"

"p14=4.875/(80.0*bf2/1000000)" /* REBURP pulse length  */
"spw3=plw2*(pow((p3*1.97/p14)/0.0798,2))"   /* REBURP power level  */
"spoffs3=0"
"spw7=plw2*(pow((p3*1.97/p14)/0.0798,2))"   /* REBURP power level  */
"spoffs7=bf2*((cnst21-cnst22)/1000000)"   /* shift from CAB to CO   */

/*******************************************************************/


"p22=p21*2"

"d0=0"
"in0=inf1/2"
"d10=0"
"in10=inf2/2"
"in20=inf2/2"

"d22=13m"
"d20=d22"


"DELTA=p16+d16+50u+d0*2"

"DELTA1=p16+d16-p1*0.78+de+8u"
"DELTA2=d4-p25*0.5"
"DELTA3=d24-p19-d16"
"DELTA4=d4-p16-d16-p25*0.5"
"DELTA5=d22-p16-d16-50u"



"acqt0=0"
baseopt_echo

aqseq 312

1 ze
  d11 pl12:f2
2 d11 do:f2 
# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  30u fq=0:f1

3 d12 
/**************************************/
/*   1H (t1) editing  & NOE mixing    */
/**************************************/
  (p1 ph8:r)
  d0
  (center (p25:sp25 ph1):f2 (p22 ph1):f3 )
  d0
  (p2 ph1)
  (center (p25:sp25 ph1):f2 (p22 ph1):f3 )
  (p1 ph9)

  p16:gp6
  d16

   d8

/**************************************/
/*   H-C transfer                     */
/**************************************/
  (p1 ph1)
  DELTA2 
  (center (p2 ph1) (p25:sp25 ph6):f2 )
  DELTA2 pl2:f2
  (p1 ph2) 
  p16:gp5
  d16

  (p3 ph3):f2

/**************************************/
/*   13C editing                      */
/**************************************/
#   ifdef CT
  d10 
  (center (p2 ph1) (p14:sp7 ph1):f2 (p22 ph1):f3 )
  DELTA5

  50u UNBLKGRAD
  p16:gp1*EA
  d16 pl2:f2
  10u
  (p14:sp3 ph4):f2
  d20
  (p14:sp7 ph1):f2
  10u pl2:f2
#   else

  d10 
  (center (p2 ph1) (p14:sp7 ph1):f2 (p22 ph1):f3 )
  d10

  50u UNBLKGRAD
  p16:gp1*EA
  d16 pl2:f2
  (p4 ph4):f2
  (p14:sp7 ph1):f2
  DELTA pl2:f2
#   endif /*CT*/
/**************************************/
/*   SE 1H-13C back transfer          */
/**************************************/
  (center (p1 ph1) (p3 ph4):f2 )
  p19:gp3
  d16
  DELTA3
  (center (p2 ph1) (p4 ph1):f2 )
  DELTA3
  p19:gp3
  d16 pl2:f2
  (center (p1 ph2) (p3 ph5):f2 )
  p16:gp4
  d16
  DELTA4 
  (center (p2 ph1) (p25:sp25 ph1):f2 )
  DELTA4
  p16:gp4
  d16
  (p1 ph1)
  DELTA1
  (p2 ph1)
  4u
  p16:gp2
  d16 pl12:f2
  4u BLKGRAD
/**************************************/
/*   Detection                        */
/**************************************/
  go=2 ph31 cpd2:f2
  d11 do:f2 mc #0 to 2 
     F1PH(calph(ph8, +90), caldel(d0, +in0))
     F2EA(calgrad(EA) & calph(ph5, +180), caldel(d10, +in10) &caldel(d20, -in20) &  calph(ph3, +180) & calph(ph6, +180) & calph(ph31, +180))
exit
   

ph1=0 
ph2=1
ph3=0 2
ph4=0 0 0 0 2 2 2 2
ph5=1 1 1 1 3 3 3 3
ph6=0
ph8=0 0 2 2
ph9=0
ph29=0
ph31=0 2 2 0 2 0 0 2
  

;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;sp3: f2 channel - shaped pulse 180 degree
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p14: f2 channel - 180 degree shaped pulse for inversion
;p16: homospoil/gradient pulse                [1 msec]
;p19: gradient pulse 2                        [500 usec]
;p22: f3 channel - 180 degree high power pulse
;p28: f1 channel - trim pulse
;d0 : incremented delay (2D)                  [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                      [30 msec]
;d16: delay for homospoil/gradient recovery
;d22: CT_delay/2 = 1/2JCC
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


;use gradient ratio:	gp 1 : gp 2 : gp 3 : gp 4
;			  80 : 20.1 :   11 :   -5    for C-13
;			  80 :  8.1 :   11 :   -5    for N-15

;for z-only gradients:
;gpz1: 80%
;gpz2: 20.1% for C-13, 8.1% for N-15
;gpz3: 11%
;gpz4: -5%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
