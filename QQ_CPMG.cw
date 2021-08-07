/*  QQ CPMG (Chris Waudby Aug 2021)
    based on 13CH3_1H_TQ_CPMG_ty_800_cpz

    This pulse sequence will allow one to perform the following experiment:

    2D 1H/13C to measure exchange via 13CH3 CPMG 

    Assumes that sample is specifically 13CH3 labeled

        1H: O1 on methyls (0.8 ppm)
            pwh = p1 1H pw90 @ power level pl1 highest power
            pwh_cp = p15 1H pw90 @ power level pl15 for CPMG pulses

       13C: O2 centre at 20 ppm 
            pwc = p2 13C pw90 @ power level pl2 highest power
            power level pl21 is used for 13C decoupling.

        2H: O4 centre at ~0.8 ppm (methyl 1H resonances)
            pwd = p4 2H pw90 @ power level pl4 high power
            power level pl41 is used for 2H waltz-16 decoupling.
            Note: 2H decoupling only necessary for measuring nuCPMG < 100 Hz

       Recommend: use -Dcomp_flg to apply composite 180o CPMG pulses

       Written by T. Yuwen on Nov 21, 2016
*/

prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

/***********************/
/*   Define pulses     */
/***********************/
define pulse dly_pg1     /* Messerle purge pulse */
  "dly_pg1=2m"
define pulse dly_pg2     /* Messerle purge pulse */
  "dly_pg2=3.4m"
define pulse pwh
  "pwh=p1"               /* 1H hard pulse at power level p1 (tpwr) */
define pulse pwc
  "pwc=p2"               /* 13C pulse at power level pl2 (dhpwr) */

#ifdef Ddec
define pulse pwd
  "pwd=p4"               /* 2H pulse at power pl4 */
#endif

define pulse pwh_cp      /* 1H CPMG pulse power level */
  "pwh_cp=p15"

/************************/
/*   Define delays      */
/************************/
"in0=inf2/2"
"d11=30m"

/************************/
/*   Define f1180       */
/************************/
  "d0=larger((in0)/2 - 2.0*pwh_cp, 2e-7)"

define delay taua
  "taua=d3"              /* d3 ~ 1.8-2ms ~ 1.0s/(4*125.3)"  ~ 1 / 4J(CH) */
define delay taub
  "taub=d4"              /* d4 = 1/4JCH exactly */ 

define delay time_T2
  "time_T2=d6"           /* CPMG duration <= 40 ms  */

/*************************************/
/* Define parameters related to CPMG */
/*************************************/
define delay tauCPMG
define delay tauCPMG1
define delay tauCPMG2

define list<loopcounter> ncyc_cp=<$VCLIST>

"cnst31=1/time_T2"       /* minimum nuCPMG for given time_T2 */

/******************************************************/
/* Assign cnsts to check validity of parameter ranges */
/******************************************************/

"cnst15=plw15"           /* tpwrcp - power level for 1H CPMG pulses */
"cnst21=plw21"           /* dpwr pl21 - set max at 2.0W */


/****************************/
/* Initialize loop counters */
/****************************/
"l1=0"
"l2=0"
"l3=0"

"acqt0=0"              /* select 'DIGIMOD = baseopt' to execute */

aqseq 312

1 ze
#ifdef Ddec
; d11 LOCKDEC_ON /* Not required for AvanceIII-HD */
  50u LOCKH_ON
  d11 H2_PULSE
  2u pl41:f4
#endif /*Ddec*/

2 d11 do:f2
/****************************************/
/* Continue to check run time variables */
/****************************************/
  "l2 = (trunc(ncyc_cp[l1] + 0.3))"

  if "l2 > 80"
  {
    2u
    print "error: ncyc_cp must be < 81"
    goto HaltAcqu
  }

  if "l2 != 0" 
  {
    "tauCPMG = time_T2*0.25/l2"
    if "tauCPMG < 124.0u"
    {
      2u
      print "error: tauCPMG < 124u too short for CPMG"
      goto HaltAcqu
    }
  }
  else 
  {
    "tauCPMG = time_T2*0.25"
  }

  "tauCPMG1 = tauCPMG - pwh_cp*2.0"
#ifdef comp_flg
  "tauCPMG2 = tauCPMG - pwh_cp*2.0"
#else
  "tauCPMG2 = tauCPMG - pwh_cp"
#endif

  20u pl1:f1 pl2:f2

#ifdef Ddec
  d11 H2_LOCK
  6m LOCKH_OFF
#endif /*Ddec*/

/******************/
/* Messerle purge */
/******************/
  20u pl11:f1
  (dly_pg1 ph26):f1
  20u
  (dly_pg2 ph27):f1
  20u pl10:f1

 
; off-resonance presat
30u fq=cnst10(bf hz):f1
30u pl10:f1
d1 cw:f1 ph26
4u do:f1
30u fq=0:f1


#ifdef Ddec
  50u LOCKH_ON
  15u H2_PULSE
#endif

/****************************************/
/* Destroy 13C equlibrium magnetization */
/****************************************/
  (pwc ph26):f2

#ifdef Ddec
  20u UNBLKGRAMP
#else
  20u UNBLKGRAD
#endif

  2u
  p50:gp0
  d16

/***********************/
/* Create TQ coherence */
/***********************/
  20u pl15:f1
  2u rpp11 rpp12 rpp13 rpp14 rpp15

  (pwh_cp ph1):f1

  2u
  p51:gp1
  d16

  "DELTA = taua - 2u - p51 - d16 - pwh_cp*2.0/PI"
  DELTA

  (center (pwh_cp*2 ph1):f1 (pwc*2 ph26):f2)

  2u
  p51:gp1
  d16

  "DELTA = taua - 2u - p51 - d16"
  DELTA

  (pwc ph3):f2

  2u
  p52:gp2
  d16

  "DELTA = taub - 2u - p52 - d16"
  DELTA

  (center (pwh_cp*2 ph1):f1 (pwc*2 ph26):f2)

  2u
  p52:gp2
  d16

#ifdef Ddec
  "DELTA = taub - 2u - p52 - d16 - 2u - pwd - 2u - 2u"
#else
  "DELTA = taub - 2u - p52 - d16"
#endif
  DELTA

#ifdef Ddec
  2u pl4:f4
  (pwd ph27):f4
  2u pl41:f4
  (2u cpds4 ph26):f4
#endif

;  (pwc ph26):f2
  (pwh_cp ph1):f1

/*********************************/
/* The first half of CPMG period */
/*********************************/
  if "l2 == 1" {
    tauCPMG1
    (pwh_cp ph29 pwh_cp*2 ph26 pwh_cp ph29):f1
    tauCPMG1
  }

  if "l2 > 1" {
3   tauCPMG2
#ifdef comp_flg
    (pwh_cp ph14 pwh_cp*2 ph11 pwh_cp ph14):f1
#else
    if "l2%2 == 0" {
      (pwh_cp*2 ph15):f1
    }
    else {
      if "abs((l2%8)-4) > 2" {
        (pwh_cp*2 ph13):f1
      }
      else {
        (pwh_cp*2 ph11):f1
      }
    }
#endif
    tauCPMG2 ipp11 ipp12 ipp13 ipp14 ipp15
    lo to 3 times l2
  }

/*****************************/
/* The central 180o 1H pulse */
/*****************************/
  (center (pwh_cp*2 ph2):f1 (pwc*2 ph25):f2 )
  
/**********************************/
/* The second half of CPMG period */
/**********************************/
  if "l2 == 1" {
    tauCPMG1
    (pwh_cp ph27 pwh_cp*2 ph26 pwh_cp ph27):f1
    tauCPMG1
  }

  if "l2 > 1" {
4   tauCPMG2 dpp11 dpp12 dpp13 dpp14 dpp15
#ifdef comp_flg
    (pwh_cp ph13 pwh_cp*2 ph12 pwh_cp ph13):f1
#else
    if "l2%2 == 0" {
      (pwh_cp*2 ph15):f1
    }
    else {
      if "abs((l2%8)-4) > 2" {
        (pwh_cp*2 ph14):f1
      }
      else {
        (pwh_cp*2 ph12):f1
      }
    }
#endif
    tauCPMG2
    lo to 4 times l2
  }

  (pwh_cp ph26):f1
;  (pwc ph3):f2

#ifdef Ddec
  2u do:f4
  2u pl4:f4
  (pwd ph29):f4
#endif

  2u
  p53:gp3
  d16

#ifdef Ddec
  "DELTA = taub - 2u - 2u - pwd - 2u - p53 - d16"
#else
  "DELTA = taub - 2u - p53 - d16"
#endif
  DELTA

/********/
/* HMQC */
/********/
  (pwc*2 ph26):f2

  d0
  (pwh_cp ph29 pwh_cp*2 ph26 pwh_cp ph29):f1
  d0

  2u
  p53:gp3
  d16

  "DELTA = taub - 2u - p53 - d16"
  DELTA

  (pwc ph4):f2

  "DELTA = pwc*2.0"
  DELTA

  (pwh_cp ph27):f1

  2u
  p54:gp4
  d16

/****************************************************************/
/* C->H back transfer, use wtg_flg for better water suppression */
/****************************************************************/
  20u pl1:f1
  
  (pwh ph26):f1

  2u
  p57:gp7
  d16
  "DELTA = taua - 2u - p57 - d16 - p10 - 1u - larger(pwh,pwc) - pwh*2.0/PI"
  DELTA
  (p10:sp10 ph28):f1
  1u pl1:f1
  (center (pwh*2 ph26):f1 (pwc*2 ph27):f2 )
  1u
  (p10:sp10 ph28):f1
  "DELTA = taua - p57 - d16 - p10 - 1u - larger(pwh,pwc) - 2*pwc - 8u"
  DELTA
  p57:gp7
  d16 

#ifdef Ddec
  4u BLKGRAMP
#else
  4u BLKGRAD
#endif

  (pwc ph26):f2
  (pwc ph5):f2

  4u pl21:f2             /* lower power for 13C decoupling */

/********************************/
/* Signal detection and looping */
/********************************/
  go=2 ph31 cpds2:f2
  d11 do:f2 mc #0 to 2
    F1QF(iu1)
    F2PH(ru1 & ip4, id0); & ip31*2)
;    F1QF(calclc(l1,1))
;    F2PH(calph(ph4,-90), caldel(d0,+in0) & calph(ph31,+180))

#ifdef Ddec
  d11 H2_LOCK
  d11 LOCKH_OFF
; d11 LOCKDEC_OFF        /* use statement for earlier hardware */
#endif

HaltAcqu, 1m
exit

ph0=1
ph1=(6) 0 1 2 3 4 5
ph2=1 3
ph3={{0}*6}^2
ph4={{0}*12}^2
ph5=0 2
ph11={0 1 0 1 1 0 1 0}^2
ph12={0 3 0 3 3 0 3 0}^2
ph13=ph12+ph0
ph14=ph11-ph0
ph15={{{0 2}^2}^2}^2
ph25=0
ph26=0
ph27=1
ph28=2
ph29=3
ph31={{{0 2}*3}^2}^2

;pl1 : tpwr - power level for pwh
;pl2 : dhpwr - power level for 13C pulse pwc (p2)
;pl4  : power level for 2H high power pulses
;pl10 : tsatpwr - power level for presat
;pl11 : tpwrmess - power level for Messerle purge
;pl15 : power level for 1H CPMG pulses pwh_cp
;pl21 : dpwr - power level for  13C decoupling cpd2
;pl41 : power level for 2H waltz decoupling
;p10 : 1000usec water flip-back
;sp10 : water flip-back (on H2O)
;spw14 : power level for eburp1 pulse
;spnam14: eburp1 pulse on water
;p1 : pwh
;p2 : pwc
;p4 : 2H high power pulse
;p14 : eburp1 pulse width, typically 7000u
;p15 : 1H pw for CPMG pulses
;p50 : gradient pulse 50                                [1000 usec]
;p51 : gradient pulse 51                                [400 usec]
;p52 : gradient pulse 52                                [600 usec]
;p53 : gradient pulse 53                                [300 usec]
;p54 : gradient pulse 54                                [500 usec]
;p55 : gradient pulse 55                                [300 usec]
;p56 : gradient pulse 56                                [500 usec]
;p57 : gradient pulse 57                                [800 usec]
;pcpd2 : 13C pulse width for 13C decoupling
;pcpd4 : 2H pulse width for 2H decoupling
;d1 : Repetition delay D1
;d3 : taua ~1/(4*JCH)  ~1.8-2ms
;d4 : taub - set to 1/4JHC = 2.0 ms
;d6 : time_T2 CPMG duration <= 40ms
;d11 : delay for disk i/o, 30ms
;d16 : gradient recovery delay, 200us
;d19: delay for binomial water suppression, = 1/(4*|cnst1|)
;cpd2 : 13C decoupling during t2 according to program defined by cpdprg2
;cpd4 : 2H decoupling during t1
;cpdprg2 : 13C decoupling during t2
;cpdprg4 : 2H decoupling during t1
;cnst1 : offset of methyls from water (0.8 ppm - 4.7 ppm, in Hz)
;cnst4 : power level for 2H high power pulses
;cnst10: water frequency for presat
;cnst14 : power level in w for 1H selective pulse
;cnst15 : power in w for 1H CPMG pulses
;cnst21 : power in w for 13C dec
;cnst31 : minimum nuCPMG for given time_T2
;cnst41 : power level for 2H decoupling
;l1 : counter for the ncyc_cp values for cpmg
;l2 : actual value of ncyc_cp
;l3 : flag to switch between IP and AP for -Dcal_flg
;inf1 : 1/SW(X) = 2*DW(X)
;in0 : 1/(2*SW(x))=DW(X)
;nd0 : 2
;ns : 6*n
;FnMODE : States-TPPI, TPPI, States

;for z-only gradients:
;gpz0: 20%
;gpz1: 25%
;gpz2: 30%
;gpz3: -25%
;gpz4: 50%
;gpz5: -40%
;gpz6: -75%
;gpz7: -80%

;use gradient files:
;gpnam0: SMSQ10.32
;gpnam1: SMSQ10.32
;gpnam2: SMSQ10.32
;gpnam3: SMSQ10.32
;gpnam4: SMSQ10.32
;gpnam5: SMSQ10.32
;gpnam6: SMSQ10.32
;gpnam7: SMSQ10.32

;zgoptns: Dfsat Dmess_flg Dfscuba Dwater_flg Dwtg_flg Df1180 DDdec Dcomp_flg Dmq_flg Dcal_flg
