/* 1H TQ relaxation measurement in methyl groups
based on 2016 Yuwen sequence

    Assumes that sample is specifically 13CH3 labeled

        1H: O1 on methyls (0.8 ppm)
            pwh = p1 1H pw90 @ power level pl1 highest power

       13C: O2 centre at 20 ppm 
            pwc = p2 13C pw90 @ power level pl2 highest power
            power level pl21 is used for 13C decoupling.

        2H: O4 centre at ~0.8 ppm (methyl 1H resonances)
            pwd = p4 2H pw90 @ power level pl4 high power
            power level pl41 is used for 2H waltz-16 decoupling.
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
  "pwc=p3"               /* 13C pulse at power level pl2 (dhpwr) */
define pulse pwd
  "pwd=p30"               /* 2H pulse at power pl4 */

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

/*************************************/
/* Define parameters related to CPMG */
/*************************************/
define delay tauEcho
  "tauEcho=666.7u"

define list<loopcounter> ncyc=<$VCLIST>


/****************************/
/* Initialize loop counters */
/****************************/
"l1=0"
"l2=0"
"l3=0"

"acqt0=0"              /* select 'DIGIMOD = baseopt' to execute */

aqseq 312

1 ze
; d11 LOCKDEC_ON /* Not required for AvanceIII-HD */
  50u LOCKH_ON
  d11 H2_PULSE
  2u pl17:f4

2 d11 do:f2
  "l2 = (trunc(ncyc[l1] + 0.3))"

  20u pl1:f1 pl2:f2

  d11 H2_LOCK
  6m LOCKH_OFF

/******************/
/* Messerle purge */
/******************/
  20u pl11:f1
  (dly_pg1 ph26):f1
  20u
  (dly_pg2 ph27):f1

 
; off-resonance presat
30u fq=cnst10(bf hz):f1
30u pl9:f1
d1 cw:f1 ph26
4u do:f1
30u fq=0:f1
20y pl1:f1

  50u LOCKH_ON
  15u H2_PULSE

/****************************************/
/* Destroy 13C equlibrium magnetization */
/****************************************/
  (pwc ph26):f2

  20u UNBLKGRAMP

  2u
  p50:gp0
  d16

/***********************/
/* Create TQ coherence */
/***********************/

  (pwh ph1):f1

  2u
  p51:gp1
  d16

  "DELTA = taua - 2u - p51 - d16 - pwh*2.0/PI"
  DELTA

  (center (pwh*2 ph1):f1 (pwc*2 ph26):f2)

  2u
  p51:gp1
  d16

  "DELTA = taua - 2u - p51 - d16"
  DELTA

  (pwc ph26):f2

  2u
  p52:gp2
  d16

  "DELTA = taub - 2u - p52 - d16"
  DELTA

  (center (pwh*2 ph1):f1 (pwc*2 ph26):f2)

  2u
  p52:gp2
  d16

  "DELTA = taub - 2u - p52 - d16 - 2u - pwd - 2u - 2u"
  DELTA

  2u pl4:f4
  (pwd ph27):f4
  2u pl17:f4
  (2u cpds4 ph26):f4

  (pwc ph26):f2

/*************/
/* Hahn echo */
/*************/
  (pwh ph1):f1
  "tauEcho=l2*666.7u - phw*2.6366"
  tauEcho
  (pwh ph29 pwh*2 ph26 pwh ph29):f1
  tauEcho

  (pwh ph26):f1
  (pwc ph3):f2

  2u do:f4
  2u pl4:f4
  (pwd ph29):f4

  2u
  p53:gp3
  d16

  "DELTA = taub - 2u - 2u - pwd - 2u - p53 - d16"
  DELTA

/********/
/* HMQC */
/********/
  (pwc*2 ph26):f2

  d0
  (pwh ph29 pwh*2 ph26 pwh ph29):f1
  d0

  2u
  p53:gp3
  d16

  "DELTA = taub - 2u - p53 - d16"
  DELTA

  (pwc ph4):f2

  "DELTA = pwc*2.0"
  DELTA

  (pwh ph27):f1

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

  4u BLKGRAMP

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
ph3={{0}*6}^2
ph4={{0}*12}^2
ph5=0 2
ph26=0
ph27=1
ph28=2
ph29=3
ph31={{{0 2}*3}^2}^2

;pl1 : tpwr - power level for pwh
;pl2 : dhpwr - power level for 13C pulse pwc (p2)
;pl4  : power level for 2H high power pulses
;pl9 : tsatpwr - power level for presat
;pl11 : tpwrmess - power level for Messerle purge
;pl15 : power level for 1H CPMG pulses pwh_cp
;pl21 : dpwr - power level for  13C decoupling cpd2
;pl17 : power level for 2H waltz decoupling
;p10 : 1000usec water flip-back
;sp10 : water flip-back (on H2O)
;spw14 : power level for eburp1 pulse
;spnam14: eburp1 pulse on water
;p1 : pwh
;p3 : pwc
;p30 : 2H high power pulse
;p14 : eburp1 pulse width, typically 7000u
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
;d11 : delay for disk i/o, 30ms
;d16 : gradient recovery delay, 200us
;cpd2 : 13C decoupling during t2 according to program defined by cpdprg2
;cpd4 : 2H decoupling during t1
;cpdprg2 : 13C decoupling during t2
;cpdprg4 : 2H decoupling during t1
;cnst10: water frequency for presat
;l1 : counter for the ncyc_cp values for cpmg
;l2 : actual value of ncyc_cp
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

