;4D HC(CC)CONH TOCSY
; for sidechain assignment
; adapted from IBS library, Chris Waudby Feb 2020
; using wavemaker (wvm)
;
; 1H(ali) [t1] --> 13C(ali) [t2] --> 15N [t3] --> 1H [t4]
;
;BB_HCCONH_TOCSY
;BB 06/12/2016

prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"
"p17=260u"
"p19=600u"


/*******************************************************************/
/*   calculation of TOCSY loop number                              */
/*******************************************************************/
; TOCSY mixing time = l1 x 2.9 ms  (p9*115.112)

/*******************************************************************/
/*   calculation of shaped 1H pulse parameters                     */
/*******************************************************************/
"p42=4.875/(cnst2*bf1/1000000)" /* REBURP pulse length  */
"spw26=plw1*(pow((p1*1.97/p42)/0.0798,2))"   /* REBURP power level  */
"spoff26=bf1*(cnst1/1000000)-o1" /* REBURP offset */
;"spoal26=0.5"

"p43=4.6/(cnst2*bf1/1000000)" /*  EBURP pulse length   */
"spw28=plw1*(pow((p1*1.04/p43)/0.06103,2))"   /* EBURP power level  */
"spoff28=bf1*(cnst1/1000000)-o1" /*  EBURP offset */
"spw29=plw1*(pow((p1*1.04/p43)/0.06103,2))"   /* EBURP_TR power level  */
"spoff29=bf1*(cnst1/1000000)-o1" /*  EBURP_TR offset */

"p44 =p1*8.0"    /* BIP pulse length  */
"spoff30=0.0" /*  BIP offset */
"spw30=plw1"  /* BIP power level  */

"p25=pcpd1"

/*******************************************************************/
/*   calculation of shaped 13C pulse parameters                    */
/*******************************************************************/
"cnst23 = (sfo2-bf2)*1000000/bf2"    /*  Caliph  frequency offset   */
"cnst22 = cnst23-(39.0-54.0)"     /*  CA frequency offset   */
"cnst21 = cnst23-(39.0-175.0)"     /*  CO frequency offset   */

"p9=25u"
"plw15=plw2*(pow((p3/p9),2))"   /* DIPSI-2 pulse power level */

/*******************************************************************/
/*   calculation of shaped 15N pulse parameters                    */
/*******************************************************************/
"p50 =500u"    /* BIP pulse length  */
"p50 = p50 +cnst11 - cnst11" ; TODO this looks odd!
"spoff50=0.0" /*  BIP offset */
"spw50=plw3*(pow((p21*8/p50),2))"   /* BIP power level  */

/*******************************************************************/
/*   DELAYS                                                        */
/*******************************************************************/
"d3=1.1m"	
"d4=1.7m"			;tau a
"d21=13.4m"			;T N
"if (l0 == 1) { d22=3.6m; } else { d22=4.4m; }"   /* CA-CO transfer delay 1/4J  */
"d24=4.4m"			;tau d
"d25=5.5m"			;tau f
"d26=2.7m"			;tau g
"d27=14m"			;tau e
"d0=3u"

"DELTA1=d26-p17-d16-p43*0.5-p42*0.5"
"DELTA2=d26-p16-d16-p42*0.5"
"DELTA3=d27-d24+4u"
"DELTA4=d4+p14+d0*2"
"TAU=d3+p2+d0*2-4u"
"DELTA5=p16+d16+8u"
"DELTA6=d22-p20*0.5"
"DELTA9=d22"
"DELTA7=d21-d26-p16-d16-p14-p44+p21*4/PI"
"DELTA8=d26-p14"

/*******************************************************************/
/*   time incremennts in 1H dimension                              */
/*******************************************************************/
"d0=3u"
"in0=inf1/2"

/*******************************************************************/
/*   time incremennts in 13C dimension                             */
/*******************************************************************/
"d20=3u"
"in20=inf2/2"

/*******************************************************************/
/*   time increments in 15N dimension                              */
/*******************************************************************/

"d10=3u"
"in10=inf3/2"
"d29=p43" 
"d30=d21+3u"

"FACTOR2=d30*10000000*2/td3"
"in30=FACTOR2/10000000"

"if ( in30 > in10 ) { in29 = 0; } else { in29=in10-in30; }"
"if ( in30 > in10 ) { in30 = in10; }"



/*******************************************************************/
/*  Start of pulse sequence                                        */
/*******************************************************************/

"acqt0=0"
baseopt_echo

1 ze
  d11 pl16:f3
2 d11 do:f3
3 d11 fq=cnst23(bf ppm):f2
  d1
  50u UNBLKGRAD
  d12 pl1:f1 pl2:f2 pl3:f3
 50u fq=cnst23(bf ppm):f2   /* 13C carrier at Caliph  */

/*******************************************************************/
/* 1H->13C transfer                                                */
/*******************************************************************/
  (p1 ph12):f1
  d0
  d4
  (p14:sp3 ph1):f2 
  d0
  (p2 ph1):f1
  DELTA4
  (p1 ph13):f1

  p16:gp6
  d16 
/*******************************************************************/
/* 13C editing and back transfer                                   */
/*******************************************************************/
 (p13:sp2 ph14):f2
  d20
  (center (p14:sp5 ph1):f2 (p22 ph1):f3 )
  d3
  (p2 ph1):f1
  d20
  (p14:sp3 ph1):f2
  TAU
  (p14:sp5 ph1):f2
  4u pl2:f2
  (p13:sp8 ph13):f2

  (p1 ph2):f1

  p16:gp6
  d16 

  10u pl15:f2

/*******************************************************************/
/* 13C TOCSY using DIPSI-2                                         */
/*******************************************************************/
						;begin DIPSI2
7 (p9*3.556 ph23):f2
  (p9*4.556 ph25):f2
  (p9*3.222 ph23):f2
  (p9*3.167 ph25):f2
  (p9*0.333 ph23):f2
  (p9*2.722 ph25):f2
  (p9*4.167 ph23):f2
  (p9*2.944 ph25):f2
  (p9*4.111 ph23):f2
  
  (p9*3.556 ph25):f2
  (p9*4.556 ph23):f2
  (p9*3.222 ph25):f2
  (p9*3.167 ph23):f2
  (p9*0.333 ph25):f2
  (p9*2.722 ph23):f2
  (p9*4.167 ph25):f2
  (p9*2.944 ph23):f2
  (p9*4.111 ph25):f2

  (p9*3.556 ph25):f2
  (p9*4.556 ph23):f2
  (p9*3.222 ph25):f2
  (p9*3.167 ph23):f2
  (p9*0.333 ph25):f2
  (p9*2.722 ph23):f2
  (p9*4.167 ph25):f2
  (p9*2.944 ph23):f2
  (p9*4.111 ph25):f2

  (p9*3.556 ph23):f2
  (p9*4.556 ph25):f2
  (p9*3.222 ph23):f2
  (p9*3.167 ph25):f2
  (p9*0.333 ph23):f2
  (p9*2.722 ph25):f2
  (p9*4.167 ph23):f2
  (p9*2.944 ph25):f2
  (p9*4.111 ph23):f2
  lo to 7 times l1
						;end DIPSI2
  d12 pl2:f2

/*******************************************************************/
/* CA->CO transfer                                                 */
/*******************************************************************/
 50u fq=cnst22(bf ppm):f2   /* 13C carrier at CA  */
 10u fq=cnst20(bf hz):f1 /* 1H carrier on water */
  d12 pl19:f1
  p25 ph1 
  d12 cpds1:f1 ph2

if "l0 ==1"   /* selective CA->CO transfer   */
{
  (p13:sp2 ph3):f2
  DELTA6
  (p20:sp10 ph1):f2
  DELTA6
  (p13:sp8 ph2):f2
}
else   /* non-selective CA-CO transfer   */
{
  
  (p13:sp2 ph3):f2
  (p14:sp4 ph1):f2
  DELTA9
  (p14:sp3 ph1):f2
  (p14:sp4 ph1):f2
  DELTA9 
  (p13:sp8 ph2):f2
}

/*******************************************************************/
/* CO->N transfer                                                 */
/*******************************************************************/
  
  p16:gp1
  d16 fq=cnst21(bf ppm):f2           /* 13C carrier at CO  */

  (p13:sp2 ph4):f2
  d24
  (p14:sp7 ph1):f2
  DELTA3
  (center (p14:sp3 ph1):f2 (p22 ph1):f3 )
  d27
  (p14:sp7 ph1):f2
  4u
  (p13:sp8 ph1):f2

  4u do:f1
  p25 ph20  /* -x */
  10u fq=0:f1 /* 1H carrier back to default */
/*******************************************************************/
/* N->CO transfer  & semi-CT 15N editing                           */
/*******************************************************************/
 (p21 ph11):f3
  (p50:sp50 ph1):f3
  d10
  (p14:sp7 ph1):f2   /* CA 180deg  */
  DELTA8   /* 1/4JNH-p14  */
  (p44:sp30 ph1)
  p16:gp2*EA
  d16
  DELTA7 
  (p14:sp3 ph1):f2   /* CO 180deg  */
  d29              /*   t2b   */
  (p50:sp50 ph1):f3
  d30  pl3:f3      /*   t2a    */

/*******************************************************************/
/* SE H-N back transfer                                            */
/*******************************************************************/
 (p43:sp28 ph1)   /* EBURP */

 (p21 ph5):f3
  p16:gp5
  d16
  DELTA2
  
  (center (p42:sp26 ph1) (p51:sp51 ph1):f3 )
  DELTA2
  p16:gp5
  d16 pl3:f3
  (p21 ph6):f3
  (p43:sp29 ph2)    /* EBURP_REV  */
/**************************************/
  p17:gp6
  d16
  DELTA1   
  (center (p42:sp26 ph2) (p51:sp51 ph2):f3 )
  DELTA1
  p17:gp6
  d16
  (p43:sp28 ph1)   /* EBURP */
/**************************************/
  DELTA5
  (p42:sp26 ph1)  /* REBURP */
  p16:gp3
  d16  pl16:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
   F1PH(calph(ph12, +90), caldel(d0, +in0)) 
   F2PH(calph(ph14, +90), caldel(d20, +in20)) 
   F3EA(calgrad(EA) & calph(ph5, +180), caldel(d10, +in10) & caldel(d29, +in29) & caldel(d30, -in30) & calph(ph11, +180) & calph(ph31, +180))


exit


ph1=0
ph2=1
ph3=0 
ph4=0 
ph5=0 0 0 0 2 2 2 2
ph6=3 3 3 3 1 1 1 1
ph7=3
ph8=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph9=2
ph10=0 0 0 0 2 2 2 2
ph11=0
ph12=0 2
ph13=1
ph14=0 0 2 2
ph20=2
ph23=0
ph25=2
ph31=0 2 2 0 2 0 0 2


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (deFLAult)
;pl2 : f2 channel - power level for pulse (deFLAult)
;pl3 : f3 channel - power level for pulse (deFLAult)
;pl16: f3 channel - power level for CPD/BB decoupling
;pl19: f1 channel - power level for CPD/BB decoupling

;                  for time reversed pulse
;p0 : f1 channel -120/60 degree high power pulse 
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p13: f2 channel -  90 degree shaped pulse
;p14: f2 channel - 180 degree shaped pulse
;p15: f2 channel - 180 degree shaped pulse (more selective for C=O) [400u @ 600MHz]
;p24: f2 channel - 180 degree shaped pulse (S/T selective) [900u @ 600MHz]
;p16: homospoil/gradient pulse                         [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p26: f1 channel -  90 degree pulse at pl19
;d0 : incremented delay (F1 in 3D)                     [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d3 : tau b : 1.1m-p14
;d4 : 1/(4J(CH)) - tau a                               [1.7 msec]
;d10: incremented delay (F2 in 3D) =  d21/2-p14/2
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d16: delay for homospoil/gradient recovery
;d21: T(N)                                             [12.4 msec]
;d23: tau c                                            [3.6 msec]
;d24: tau d                                            [4.4 msec]
;d25: tau f                                            [5.5 msec]
;d26: 1/(4J(NH)) - tau g                               [2.3 msec]
;d27: tau e                                            [12.4 msec]
;d29: incremented delay (F2 in 3D) = d21/2-p14/2-p26-d25-4u
;d30: decremented delay (F2 in 3D) = d21/2-p14/2
;cnst1: H(N) excitation frequency (in ppm)
;cnst2: H(N) excitation band width (in ppm)
;cnst9: Ser/Thr CB chemical shift (offset, in ppm) [72 ppm]
;cnst11: 15N decoupling bandwidth in detection [40 ppm]
;cnst20: water frequency (Hz)
;cnst21: CO chemical shift (offset, in ppm) [180]
;cnst22: Calpha chemical shift (offset, in ppm) [54]
;cnst23: Caliphatic chemical shift (offset, in ppm) [39]
;cnst24: CO chemical shift for CBCG discrimination (offset, in ppm) [190]
;o2p: Caliphatic chemical shift (cnst23)
;inf3: 1/SW(N) = 2 * DW(N)
;in10: 1/(4 * SW(N)) = (1/2) DW(N)
;nd10: 4
;in29: = in10
;in30: = in10
;NS: 8 * n
;DS: >= 16
;td1: number of experiments in F1 (1H)      
;td2: number of experiments in F2 (13C)
;td3: number of experiments in F3 (15N)      td2 max = 2 * d30 / in30
;FnMODE: States-TPPI (or TPPI) in F1
;FnMODE: States-TPPI (or TPPI) in F2
;FnMODE: echo-antiecho in F3
;cpds1: decoupling according to sequence defined by cpdprg1
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd1: f1 channel - 90 degree pulse for decoupling sequence
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    	gp 1 : gp 2 : gp 3 : 
;				  30 :   80 :  8.1

;for z-only gradients
;gpz1: 30%
;gpz2: 80%
;gpz3: 8.1%
;gpz4: 5%
;gpz5: -2%
;gpz7: 50%

;use gradient files:
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam5: SINE.100
;gpnam7: SINE.20


;;         WAVEMAKER -> execute: wvm_p.py
;;******************************************
;; 1H shaped pulses
;;******************************************

;cpds1(pl19):wvm:dipsi2_12p_H:dipsi2(12 ppm)

;;******************************************
;; N15 shaped pulses
;;******************************************

;sp51:wvm:reburp(40 ppm) np=1000

;cpd3(pl16):wvm:waltz16(cnst11 ppm) p90=p21 

;;******************************************
;; C13 shaped pulses
;;******************************************

;sp2:wvm:eburp2(cnst3 ppm) np=1000

;sp3:wvm:Q3(cnst3 ppm) np=1000

;sp4:wvm:Q3(cnst3 ppm, 175 ppm) np=1000 ofs=54 ppm

;sp5:wvm:Q3(cnst3 ppm, 175 ppm) np=1000 ofs=39 ppm

;sp7:wvm:Q3(cnst3 ppm, 54 ppm) np=1000 ofs=175 ppm

;sp8:wvm:eburp2_fb(cnst3 ppm) np=1000

;sp10:wvm:Q3(24 ppm, 54 ppm) Q3(24 ppm, 175 ppm) np=1000 ofs=54 ppm BS=1

