;IBS_Best_HNCA
;BH-HNCA
;using broadband 15N pulses
;with option for NUWS in 13C dimension
;ZGOPTION -DINEPT
;BB 10/01/14
;
;$CLASS=BB-Assign
;$DIM=3D
;$TYPE=BEST-HSQC
;$SUBTYPE=
;$COMMENT=

;NUWS: the vclist is applied to the 13C dimension
;Note that it does NOT require 1 to be added to the first point (like some 2D experiments)

prosol relations=<IBS>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

;define delay DELTA9
;define delay DELTA10

/*******************************************************************/
/*   calculation of shaped 1H pulse parameters                     */
/*******************************************************************/
"p41=7.2/(cnst2*bf1/1000000)" /*  PC9  pulse length  */
"spw25=plw1*(pow((p1*1.01/p41)/0.125,2))" /* PC9  power level  */
"spoff25=bf1*(cnst1/1000000)-o1"  /*  PC9  offset */
;"spoal25=0.5"

"p42=4.875/(cnst2*bf1/1000000)" /* REBURP pulse length  */
"spw26=plw1*(pow((p1*1.97/p42)/0.0798,2))"   /* REBURP power level  */
"spoff26=bf1*(cnst1/1000000)-o1" /* REBURP offset */
;"spoal26=0.5"

"p43=4.6/(cnst2*bf1/1000000)" /*  EBURP pulse length   */
"spw28=plw1*(pow((p1*1.04/p43)/0.06103,2))"   /* EBURP power level  */
"spoff28=bf1*(cnst1/1000000)-o1" /*  EBURP offset */
#   ifdef INEPT
#   else    
"spw29=plw1*(pow((p1*1.04/p43)/0.06103,2))"   /* EBURP_TR power level  */
#   endif /*INEPT*/

;"spoal28=0"
"spoff29=bf1*(cnst1/1000000)-o1" /*  EBURP_TR offset */
;"spoal29=1.0"

"p44 =p1*8.0"    /* BIP pulse length  */
"spoff30=0.0" /*  BIP offset */
"spw30=plw1"  /* BIP power level  */


/*******************************************************************/
/*   calculation of shaped 13C pulse parameters                    */
/*******************************************************************/

"p13=4.6/(60.0*bf2/1000000)" /*  EBURP pulse length   */
"spw2=plw2*(pow((p3*1.04/p13)/0.06103,2))"   /* EBURP power level  */
"spw8=plw2*(pow((p3*1.04/p13)/0.06103,2))"   /* EBURP_TR power level  */
;"spoal2=0"
;"spoal8=1"

"p14=4.875/(70.0*bf2/1000000)" /* REBURP pulse length  */
"spw3=plw2*(pow((p3*1.97/p14)/0.0798,2))"   /* REBURP power level  */
;"spw5=plw2*(pow((p3*1.97/p14)/0.0798,2))"   /* REBURP power level  */
"spw7=plw2*(pow((p3*1.97/p14)/0.0798,2))"   /* REBURP power level  */


;"cnst21 = 173.0"    /*  CO  frequency offset   */
;"cnst22 = 54.0"     /*  CA frequency offset   */

"spoff2=0"
"spoff3=0"
"spoff5=bf2*((cnst22-cnst21)/1000000)" /*  shift from CO to CA   */
"spoff7=bf2*((cnst21-cnst22)/1000000)"   /* shift from CA to CO   */
"spoff8=0"
"spoff9=0"

/*******************************************************************/
/*   calculation of shaped 15N pulse parameters                    */
/*******************************************************************/
"p50 =500u"    /* BIP pulse length  */
"spoff50=0.0" /*  BIP offset */
"spw50=plw3*(pow((p21*8/p50),2))"   /* BIP power level  */

"p51=4.875/(40*bf3/1000000)" /* REBURP pulse length  */
"spw51=plw3*(pow((p21*1.97/p51)/0.0798,2))"   /* REBURP power level  */
"spoff51=0.0"                /* REBURP offset */
;"spoal51=0.5"

/*******************************************************************/
/*   DELAYS                                                        */
/*******************************************************************/
"d11=30m"

"d25=2.5m"
"d26=2.7m"

#   ifdef INEPT
"d23=12m"
"DELTA4=d23-d26-p14-p44+p21*4/PI"
#   else    
"d23=12m-p43"
"DELTA4=d23-d26-p16-d16-p14-p44+p21*4/PI"
#   endif /*INEPT*/
"d24=12m"

"DELTA1=d25-p41*0.5-p42*0.5+p50*0.5"
"DELTA2=d26-p16-d16-p42*0.5"
"DELTA3=d26-p17-d16-p43*0.5-p42*0.5"
"if (p50 > p14) DELTA5=p50-p14"
"if (p50 <= p14) DELTA5=0"
"DELTA6=d24-p51*0.5"
"DELTA7=d24-p51*0.5-d26-p44"
"DELTA8=p16+d16+de+8u"
"DELTA9=d26-p14"



/*******************************************************************/
/*   time incremennts in 13CO dimension                            */
/*******************************************************************/
"d0=3u"
"in0=inf1/2"

/*******************************************************************/
/*   time incremennts in 15N dimension                             */
/*******************************************************************/

;"d10=10u"
;"in10=inf2/2"
;"cnst30=0.5*(td2*0.5-1)*inf2/1000000"
;"if(cnst30 > d23) in29=(cnst30-d23)/(td2*0.5-1)"
;"if(cnst30 <= d23) in29=0"
;"d30=d23+10u"   /*  t2a  */
;"if(cnst30 > d23) in30=in10-in29"
;"if(cnst30 <= d23) in30=in10"

"d10=3u"
"d30=d23+3u"
#   ifdef INEPT
"d29=0"   /*  t2b  */
#   else
"d29=p43"   /*  t2b  */
#   endif /*INEPT*/
"in10=inf2/2"

"FACTOR2=d30*10000000*2/td2"
"in30=FACTOR2/10000000"

"if ( in30 > in10 ) { in29 = 0; } else { in29=in10-in30; }"
"if ( in30 > in10 ) { in30 = in10; }"

/*******************************************************************/

#ifdef NUWS
define loopcounter dsFlag
define list<loopcounter> nuwslist=<$VCLIST>
"dsFlag=1"
#endif /* NUWS */


"acqt0=0"
baseopt_echo

aqseq 312

1 d11 ze
2 d11 do:f3
3 5u
  d1 pl1:f1 pl2:f2 pl3:f3
  50u UNBLKGRAD
/**************************************/
/*   H-N transfer                     */
/**************************************/
  (p41:sp25 ph1)    /*  PC9  */
 
  DELTA1
  (center (p42:sp26 ph2) (p50:sp50 ph1):f3 )

  DELTA1
  (p41:sp25 ph2):f1   /*  PC9  */

  p16:gp3
  d16 pl3:f3
/**************************************/
/*   15N-13CA transfer                */
/**************************************/
 30u fq=cnst22(bf ppm):f2  /* F2 carrier at CA  */

  (p21 ph1):f3
  DELTA6
  (center (p14:sp3 ph1):f2 (p51:sp51 ph1):f3 )  /* CO ,N 180deg  */
  DELTA7 pl3:f3
  (p44:sp30 ph1)
  d26
  (p21 ph1):f3

  p16:gp4
  d16
/************************************************/
/*  13CA editing                               **/
/************************************************/
; total transverse time for evolution of 1JCC
; larger(p44,p14,p50)*2 + p14 + 20u [+ p13]?
  (p13:sp2 ph14):f2   /* CA 90deg Exc  */
  d0
  (center (p44:sp30 ph1):f1 (p14:sp7 ph1):f2 (p50:sp50 ph1):f3 )   /* H, CO ,N 180deg  */
  d0
  (p14:sp3 ph1):f2
  4u
  (center (p44:sp30 ph1):f1 (p14:sp7 ph1):f2 (p50:sp50 ph1):f3 )   /* H, CO ,N 180deg  */
  10u
  (p13:sp8 ph1):f2   /* CA 90deg FB  */
  
 
  p17:gp3
  d16 pl3:f3
/************************************************/
/*  15N-13CA INEPT  with 15N editing           **/
/************************************************/
#   ifdef INEPT
  (p21 ph10):f3
  (p50:sp50 ph1):f3
  d10
  (p14:sp7 ph1):f2   /* CO 180deg  */
  DELTA9   /* 1/4JNH-p14  */
  (p44:sp30 ph1)
  DELTA4 
  (p14:sp3 ph1):f2   /* CO 180deg  */
  d29              /*   t2b   */
  (p50:sp50 ph1):f3
  d30  pl3:f3      /*   t2a    */

#   else
  (p21 ph10):f3
  (p50:sp50 ph1):f3
  d10
  (p14:sp7 ph1):f2   /* CA 180deg  */
  DELTA9   /* 1/4JNH-p14  */
  (p44:sp30 ph1)
  DELTA4 
  p16:gp1*EA
  d16
  (p14:sp3 ph1):f2   /* CO 180deg  */
  d29              /*   t2b   */
  (p50:sp50 ph1):f3
  d30  pl3:f3      /*   t2a    */
 (p43:sp28 ph1)   /* EBURP */
#   endif /*INEPT*/

/**************************************/
/*   H-N back transfer                */
/*   INEPT version                    */
/**************************************/
#   ifdef INEPT
  (p21 ph4):f3
  p17:gp5
  d16

  (p43:sp28 ph1)   /* EBURP */
  p16:gp6
  d16
  DELTA2
  (center (p42:sp26 ph2) (p51:sp51 ph1):f3 )
  DELTA2 pl16:f3
  p16:gp6
  d16 BLKGRAD
#   else
/**************************************/
/*   H-N back transfer                */
/*   SE version                       */
/**************************************/
 (p21 ph4):f3
  p16:gp5
  d16
  DELTA2

  (center (p42:sp26 ph1) (p51:sp51 ph1):f3 )
  DELTA2
  p16:gp5
  d16 pl3:f3
  (p21 ph5):f3
  (p43:sp29 ph2)    /* EBURP_REV  */
/**************************************/
  p17:gp6
  d16
  DELTA3
  (center (p42:sp26 ph2) (p51:sp51 ph2):f3 )
  DELTA3
  p17:gp6
  d16
  (p43:sp28 ph1)   /* EBURP */
/**************************************/
  DELTA8
  (p42:sp26 ph1)  /* REBURP */
  p16:gp2
  d16  pl16:f3
  4u BLKGRAD
#   endif /*INEPT*/

/**************************************/
/*   Signal detection & looping       */
/**************************************/
  go=2 ph31 cpds3:f3

#ifdef NUWS
  if "dsFlag==0" goto 10
  zd
  "dsFlag=0"
  goto 2  ; repeat following ds (without counting it as part of vclist)
10 4u
  ; repeat acquisition block according to schedule in vclist
  lo to 2 times nuwslist
#endif /* NUWS */



  d11 do:f3 mc #0 to 2 
#ifdef NUWS
     F1PH(calph(ph14, +90) & calclist(nuwslist,1), caldel(d0, +in0)& calph(ph14, +180)& calph(ph31, +180))
#else
     F1PH(calph(ph14, +90), caldel(d0, +in0)& calph(ph14, +180)& calph(ph31, +180))
#endif
#   ifdef INEPT
  F2PH(calph(ph10, +90), caldel(d10, +in10) & caldel(d29, +in29) & caldel(d30, -in30) & calph(ph10, +180) & calph(ph31, +180))
#   else
  F2EA(calgrad(EA) & calph(ph5, +180), caldel(d10, +in10) & caldel(d29, +in29) & caldel(d30, -in30) & calph(ph10, +180) & calph(ph31, +180))
#   endif /*INEPT*/

exit


ph1=0
ph2=3 
ph3=2
ph4=0
ph5=1
ph10=0 0 2 2
ph14=0 2
ph18=0 0 0 0 2 2 2 2
ph31=0 2 2 0


;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;spnam2: Eburp2.1000
;spnam3: Reburp.1000
;spnam7: Reburp.1000
;spnam8: Eburp2tr.1000
;spoal2: 0
;spoal3: 0.5
;spoal7: 0.5
;spoal8: 1
;spnam25: Pc9_4_90.1000
;spnam26: Reburp.1000
;spnam28: Eburp2.1000
;spnam29: Eburp2tr.1000
;spnam30: Bip720,50,20.1
;spoal25: 0.5
;spoal26: 0.5
;spoal28: 0
;spoal30: 0.5
;spnam50: Bip720,50,20.1
;spnam51: Reburp.1000
;spoal50: 0.5
;spoal51: 0.5
;p1: 90 degree hard pulse (1H)
;p3: 90 degree hard pulse (13C)
;p21: 90 degree hard pulse (15N)
;p16: homospoil/gradient pulse                         [500 usec]
;p17: homospoil/gradient pulse                         [300 usec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p17: gradient pulse 3                                 [300 usec]
;p41: PC9
;p42: REBURP
;p43: EBURP
;d0 : incremented delay (F1)                           [3 usec]
;d1 : relaxation delay [200 ms]
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery [200 usec]
;d25: 1/(4J(NH)                                     
;d26: 1/(4J(NH)     
;d27: 1/(4J(NH)                            
;cnst1: H(N) excitation frequency (in ppm) [8.2 ppm]
;cnst2: H(N) excitation band width (in ppm) [3.87 ppm]
;cnst21: CO chemical shift offset (in ppm) [173 ppm]
;cnst22: Calpha chemical shift offset (in ppm) [54 ppm]
;cnst26: Call chemical shift (offset, in ppm)          [101 ppm]
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(N)) = DW(N)
;nd0: 2
;ns: 4 * n
;ds: >= 16
;td1: number of experiments in F1
;o2p: 54 ppm
;o3p: 118 ppm

;for z-only gradients:

;use gradient files:   
;gpnam1: SMSQ10.32
;gpnam2: SMSQ10.32
;gpnam3: SMSQ10.32
;gpnam4: SMSQ10.32
;gpnam5: SMSQ10.32
;gpnam6: SMSQ10.32
;gpnam7: SMSQ10.32
;gpz3: -15%
;gpz4: 30%
;gpz5: 22%
;gpz6: 60%

;zgoptns: use -DINEPT


;Processing

;PHC0(F1): 45.0



;$Id: b_trosyetf3gpsi,v 1.2.2.1.4.1 2012/01/31 17:56:18 ber Exp $
