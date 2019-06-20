/* 13CH3_1H_SQ_CPMG_lek_800_cp

This pulse sequence will allow one to perform the following experiment:

2D 1H/13C to measure exchange using 1H SQ magnetization from methyl groups

      (tau 180 tau)ncyc C 180y C (tau 180 tau)ncyc

where C refers to compensation 1H 180o pulses that compensate for the fact that
  starting from SQ coherence, different coherences are created during the
  evolution of the pulse

Assumes that sample is specifically 13CH3 labeled

    1H: O1 on methyl groups (~1.0ppm)
        pwh = p1 1H pw90 @ power level pl1 highest power
        pwh_cp = p15 1H pw90 @ power level pl15 for CPMG pulses

   13C: O2 centre at 20 ppm
        pwc = p2 13C pw90 @ power level pl2 highest power
        power level pl21 is used for 13C decoupling.

Pulse sequence has the option to use regular 180o 1H pulses or 90x240y90x
  (-Dcomp180_flg); the composite pulses are recommended

Pulse sequence has the option to begin the CPMG with equal amounts of inphase
  and antiphase (-Dipap_flg) so as to minimize the effects of different
  relaxation between the two that results from Cz. Recommend to use it - there
  is no penalty in terms of extra delays

Recommend: use -Dipap_flg -Dcomp180_flg -Dwater_flg -Df1180

This method compensates so that the number of 1H 180os is fixed. Does not
  include any 180o in the reference plane. The alternative is to set -Dref_flg
  that then includes the full number of 1H 180 in the reference plane (gives
  lower R2,eff) - not recommended

Sequence has option for reburp flag in the center of the CPMG period - not used.

The sequence uses a fixed time_T2 that is independent of the number of 1H 180o
  pulses

Sequence uses xy4 based phase cycle as simulations so that this is preferred
  over xy16

Use ncyc_max = 4*k
*/

prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

/***********************/
/* Define phases */
/***********************/
#define zero ph=0.0
#define one ph=90.0
#define two ph=180.0
#define three ph=270.0

/***********************/
/* Define pulses */
/***********************/
define pulse dly_pg1 /* Messerle purge pulse */
  "dly_pg1=5m"
define pulse dly_pg2 /* Messerle purge pulse */
  "dly_pg2=dly_pg1/1.62"
define pulse pwh
  "pwh=p1" /* 1H hard pulse at power level p1 (tpwr) */
define pulse pwc
  "pwc=p2" /* 13C pulse at power level pl2 (dhpwr) */
define pulse pwh_cp /* 1H CPMG pulse power level */
  "pwh_cp=p15"

#ifdef water_flg
  define pulse pw_sl1
    "pw_sl1=p14" /* Eburp1 pulse, ~7000 us */
#endif

#ifdef reb_flg
  define pulse pwh_reb
  "pwh_reb=4.875/(cnst8*bf1/1e6)" /* REBURP pulse length */
  "spw8=plw15*(pow((p15*2.0/pwh_reb)/0.07981,2))" /* REBURP power level */
#endif /*reb_flg*/

/************************/
/* Define delays */
/************************/
define delay hscuba /* length of 1/2 scuba delay */
  "hscuba=30m"
define delay taua
  "taua=d3" /* d3 = 1/4JHC exactly */
define delay time_T2
  "time_T2=d6" /* CPMG duration <= 40 ms */

"in0=inf1/2"
"d11=30m"
"TAU2=0.2u"

/************************/
/* Define f1180 */
/************************/
#ifdef f1180
  "d0=(in0/2)"
#else
  "d0=(0.2u/2)"
#endif

/*************************************/
/* Define parameters related to CPMG */
/*************************************/
define delay tauCPMG
define delay tauCPMG1
define list<loopcounter> ncyc_cp=<$VCLIST>

/******************************************************/
/* Assign cnsts to check validity of parameter ranges */
/******************************************************/
#ifdef fsat
  "cnst10=plw10" /* tsatpwr pl10 - set max at 0.00005W */
#endif

#ifdef mess_flg
  "cnst11=plw11" /* tpwrmess pl11 - set max at 1.0W */
#endif

#ifdef water_flg
  "cnst14=spw14" /* power level for eburp1 pulse preeceding start of sequenc e */
#endif

"cnst15=plw15" /* tpwrcp - power level for 1H CPMG pulses */
"cnst21=plw21" /* dpwr pl21 - set max at 2.0W */

/**********************/
/* Define CPMG pulses */
/**********************/
#ifdef comp180_flg
#define cpmg_11 (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1
#define cpmg_13 (pwh_cp ph14 pwh_cp*2.66667 ph13 pwh_cp ph14):f1
#define cpmg_21 (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1
#define cpmg_23 (pwh_cp ph24 pwh_cp*2.66667 ph23 pwh_cp ph24):f1
"cnst51=2.333335"
#else
#define cpmg_11 (pwh_cp*2.0 ph11):f1
#define cpmg_13 (pwh_cp*2.0 ph13):f1
#define cpmg_21 (pwh_cp*2.0 ph21):f1
#define cpmg_23 (pwh_cp*2.0 ph23):f1
"cnst51=1.0"
#endif

#define cpmg_F if "(nsdone+2)%8 < 4" {\n cpmg_11 \n}\n else {\n cpmg_13 \n}
#define cpmg_R if "(nsdone+2)%8 < 4" {\n cpmg_21 \n}\n else {\n cpmg_23 \n}

#ifndef no_compensate
  define loopcounter ncyc_max /* max value of ncyc used */
    "ncyc_max=l8"
  "DELTA8 = pwh_cp*2.0*cnst51*ncyc_max*2.0"
#endif

/************************/
/* Initialize variables */
/************************/
"l1=0"
"l2=0"
"l3=0"
"spoal8=0.5"
"spoff8=0"
"spoal14=1"
"spoff14=0"

aqseq 321

"acqt0=0"
baseopt_echo

1 ze
/******************************************************************/
/* Check validity of parameters and assign values to some of them */
/******************************************************************/
#ifdef fsat
  if "cnst10 > 0.00005" {
    2u
    print "error: tpwrmess pl10 too large; < 0.00005W !!!"
    goto HaltAcqu
  }
#endif
#ifdef mess_flg
  if "cnst11 > 1" {
    2u
    print "error: tpwrmess pl11 too large; < 1W !!!"
    goto HaltAcqu
  }
#endif
#ifdef water_flg
  if "cnst14 > 0.01" {
    2u
    print "error: power level for eburp1 pulse is too large; < 0.01W !!!"
    goto HaltAcqu
  }
#endif
if "cnst15 > 15" {
  2u
  print "error: 1H CPMG power pl15 too large; < 15W !!!"
  goto HaltAcqu
}
if "time_T2 > 40.1m" {
  2u
  print "error: time_T2 too long; < 41ms !!!"
  goto HaltAcqu
}
#ifndef no_compensate
  if "ncyc_max > 80" {
    2u
    print "error: ncyc_max too large; < 80 !!!"
    goto HaltAcqu
  }
  if "DELTA8 > 11m" {
    2u
    print "error: CPMG pulse duration too long; < 10ms !!!"
    goto HaltAcqu
  }
#endif
if "cnst21 > 6.0" {
  2u
  print "error: dpwr pl21 too large; < 2.0W !!!"
  goto HaltAcqu
}
if "aq > 64m" {
  2u
  print "error: aq is too long; < 64ms !!!"
  goto HaltAcqu
}

2 d11 do:f2
/************************/
/* Update list pointers */
/************************/
2u
"ncyc_cp.idx=l1"
2u rpp11 rpp12 rpp13 rpp14 rpp21 rpp22 rpp23 rpp24

/****************************************/
/* Continue to check run time variables */
/****************************************/
"l2 = (trunc(ncyc_cp + 0.3))"
2u
#ifdef no_compensate
  "l3 = 0"
#else /*no_compensate*/
#ifdef ref_flg
  "l3 = (trunc(ncyc_max - l2 + 0.3))"
#else /*ref_flg*/
  if "l2 > 0 " {
    "l3 = (trunc(ncyc_max - l2 + 0.3))"
  }
  else {
    "l3 = 0"
  }
#endif /*ref_flg*/
#endif /*no_compensate*/
if "l2 > 80" {
  2u
  print "error: ncyc_cp must be < 81 !!!"
  goto HaltAcqu
}
if "l2 > 0" {
  "tauCPMG = (time_T2*0.25)/l2"
#ifdef no_compensate
  "tauCPMG1 = tauCPMG - pwh_cp*cnst51"
#else
  "tauCPMG1 = tauCPMG - (DELTA8*0.25 + 0.2u*l3)/l2"
#endif
}
else {
  "tauCPMG = time_T2*0.25"
  "tauCPMG1 = 2u"
}

/**********************************/
/* 1H heating compensation period */
/**********************************/
4u pl15:f1

#if defined(ref_flg) || defined(no_compensate)
  "DELTA = 20u"
#else
  if "l2 == 0" {
    "DELTA = DELTA8 + 20u"
  }
  else {
    "DELTA = 20u"
  }
#endif

DELTA cw:f1 ph26
2u do:f1

/*************************************************/
/* Destroy residual 1H magnetization prior to d1 */
/*************************************************/
20u UNBLKGRAD

10u fq=cnst1(sfo hz):f1 /* 1H SFO1 @ tof(water) */
4u pl1:f1 /* power pl1 for 1H pulses */
(pwh ph26):f1

2u
p50:gp0*0.5
d16

(pwh ph27):f1

2u
p50:gp0
d16

4u BLKGRAD

/******************/
/* Messerle purge */
/******************/
#ifdef mess_flg
  4u pl11:f1
  (dly_pg1 ph26):f1
  2u
  (dly_pg2 ph27):f1
#endif

/*****************/
/* Presaturation */
/*****************/
#ifdef fsat
  4u pl10:f1
  d1 cw:f1 ph26
  2u do:f1
  4u pl1:f1
#ifdef fscuba
  hscuba
  (pwh ph26 pwh*2.0 ph27 pwh ph26):f1
  hscuba
#endif /*fscuba*/
#else /*fsat*/
  4u pl1:f1
  d1
#endif /*fsat*/
20u UNBLKGRAD

/**************************/
/* Water selective Eburp1 */
/**************************/
#ifdef water_flg
  2u
  (pw_sl1:sp14 ph26):f1
  2u
  2u
  p50:gp0
  d16
#endif

/****************************************/
/* Destroy 13C equlibrium magnetization */
/****************************************/
4u pl2:f2
(pwc ph26):f2

2u
p50:gp0
d16

/**************************/
/* This is the real start */
/**************************/
10u fq=0(sfo hz):f1
4u pl1:f1

(pwh ph26):f1

2u
p51:gp1
d16

"DELTA = taua*2.0 - 2u - p51 - d16 - pwh*2.0/PI"
DELTA

/*******************/
/* t1 period */
/*******************/
(pwc ph1):f2

"TAU1=larger(d0-pwh*2.0-pwc*2.0/PI,TAU2)"
TAU1

(pwh ph26 pwh*2.0 ph27 pwh ph26):f1

TAU1

(pwc ph26):f2

/************************************************/
/* Option to create 50%/50% IP/AP prior to CPMG */
/************************************************/
#ifdef ipap_flg
  if "nsdone%4 < 2" {
    "DELTA = taua + 4u"
    DELTA
    (pwh ph26 pwh*2.0 ph27 pwh ph26):f1
  }
  else {
    "DELTA = taua - 4u"
    DELTA
  }
#else
  4u
  (pwh ph26 pwh*2.0 ph27 pwh ph26):f1
#endif
4u pl15:f1

/*********************************/
/* The first half of CPMG period */
/*********************************/
  if "l2 > 0" {
3   tauCPMG1
    cpmg_F
    tauCPMG1 ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 3 times l2

  }

  if "l3 > 0" {
4   0.2u
    cpmg_F
    0.2u ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 4 times l3

}

/*****************************/
/* The central 180o 1H pulse */
/*****************************/
#ifdef reb_flg
  4u
  (pwh_reb:sp8 ph2):f1
  4u pl15:f1
#else /*reb_flg*/
  (pwh_cp*2.0 ph2):f1
#endif /*reb_flg*/

/**********************************/
/* The second half of CPMG period */
/**********************************/
  if "l3 > 0" {
5   0.2u dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    cpmg_R
    0.2u
    lo to 5 times l3

  }

  if "l2 > 0" {
6   tauCPMG1 dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    cpmg_R
    tauCPMG1
    lo to 6 times l2

  }

/**********************/
/* C->H back transfer */
/**********************/
#ifdef ipap_flg
  if "nsdone%4 >= 2" {
    4u pl1:f1
    (pwh ph26 pwh*2.0 ph29 pwh ph26):f1
    4u
  }
#endif
#ifdef ipap_flg
  "DELTA = taua - 2u - p51 - d16 - 4u - 4u - de"
#else
  "DELTA = taua*2.0 - 2u - p51 - d16 - 4u - 4u - de"
#endif
DELTA

2u
p51:gp1
d16

/************************************************************/
/* Option for 3-9-19 Watergate for better water suppression */
/************************************************************/
#ifdef wgate_flg
  4u pl1:f1
  "DELTA = 4u + 4u + de + 4u + de"
  DELTA

  2u
  p52:gp2
  d17

  (pwh*0.231 ph27):f1
  d19*2.0
  (pwh*0.692 ph29):f1
  d19*2.0
  (pwh*1.462 ph27):f1
  d19*2.0
  (pwh*1.462 ph27):f1
  d19*2.0
  (pwh*0.692 ph29):f1
  d19*2.0
  (pwh*0.231 ph27):f1

  2u
  p52:gp2
  d17
#endif

4u BLKGRAD
4u pl21:f2 /* lower power for 13C decoupling */

/********************************/
/* Signal detection and looping */
/********************************/
go=2 ph31 cpds2:f2
  d11 do:f2 mc #0 to 2
  F2QF(calclc(l1, 1))
  F1PH(calph(ph1, +90), caldel(d0, +in0) & calph(ph1, +180) & calph(ph31, +180))

HaltAcqu, 1m
exit

ph0=1
ph1=0 2
ph2=1 1 3 3 1 1 3 3 3 3 1 1 3 3 1 1
ph11=0 1 0 1
ph12=ph11-ph0
ph13=1 0 1 0
ph14=ph13-ph0
ph21=0 3 0 3
ph22=ph21+ph0
ph23=3 0 3 0
ph24=ph23+ph0
ph26=0
ph27=1
ph28=2
ph29=3
ph31=0 2

;d1: Repetition delay D1
;d3: taua - set to 1/4JHC = 2.0 ms
;d6: time_T2 CPMG duration <= 40ms
;d11: delay for disk i/o, 30ms
;d16: gradient recovery delay, 200us
;d17: gradient recovery delay for 3-9-19 watergate, 200us
;d19: delay for binomial water suppression, ~1/(4*|cnst1|)
;pl1: tpwr - power level for pwh
;pl2: dhpwr - power level for 13C pulse pwc (p2)
;pl10: tsatpwr - power level for presat
;pl11: tpwrmess - power level for Messerle purge
;pl15: power level for 1H CPMG pulses pwh_cp
;pl21: dpwr - power level for 13C decoupling cpd2
;sp14: power level for eburp1 pulse
;spnam8: Reburp.1000
;spnam14: eburp1 pulse on water
;p1: pwh
;p2: pwc
;p14: eburp1 pulse width, typically 7000us
;p15: 1H pw for CPMG pulses
;p50: gradient pulse 50 [1000 usec]
;p51: gradient pulse 51 [500 usec]
;p52: gradient pulse 52 [800 usec]
;cpdprg2: 13C decoupling program during t2 [waltz16]
;pcpd2: 13C pulse width for 13C decoupling
;cnst1: offset of water from methyls (Hz)
;cnst8: methyl H excitation bandwidth (ppm)
;vclist: variable counter list for ncyc_cp
;l8: ncyc_max (MUST BE SET PROPERLY!)
;delta8: total duration of 1H CPMG pulses
;inf1: 1/SW(X) = 2*DW(X)
;in0: 1/(2*SW(x))=DW(X)
;nd0: 2
;ns: 4*n
;FnMODE: States in F1
;FnMODE: QF in F2

;for z-only gradients:
;gpz0: 20%
;gpz1: 30%
;gpy2: 0% (Z-gradient) or 80% (XYZ-gradient)
;gpz2: 80% (Z-gradient) or 0% (XYZ-gradient)

;use gradient files:
;gpnam0: SMSQ10.32
;gpnam1: SMSQ10.32
;gpnam2: SMSQ10.32

;zgoptns: Dfsat, Dmess_flg, Dfscuba, Dwater_flg, Dwgate_flg, Df1180, Dcomp180_flg, D ipap_flg, Dreb_flg, Dref_flg, Dno_compensate
