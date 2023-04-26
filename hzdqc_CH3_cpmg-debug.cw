;H(Z/D)QC CPMG
; set d20 = relaxation time
; set td = 4x desired number of (real) points
;          (2x multiplet suppression, 2x ZQ/DQ selection)
;TODO check which is Z and D!
;
;Option for off-resonance presat (-DOFFRES_PRESAT)

prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"d11=30m"
"d12=20u"
"d13=4u"
"d21=1s/(cnst4*2)"
"d22=1s/(cnst4*8)"

"p2=p1*2"
"p4=p3*2"

"in0=inf2"
"d0=in0/2-p3*4/3.1415"


define delay tauCPMG
define delay tauCPMG1
define list<loopcounter> ncyc_cp=<$VCLIST>

define pulse pwh
  "pwh=p1" /* 1H hard pulse at power level p1 */
define pulse pwc
  "pwc=p3" /* 13C pulse at power level pl2 */
define pulse pwh_cp /* 1H CPMG pulse power level */
  "pwh_cp=p15"

define delay time_T2
  "time_T2=d20" /* CPMG duration <= 40 ms */
"TAU2=0.2u"


/**********************/
/* Define CPMG pulses */
/**********************/
;#ifdef comp180_flg
  "cnst51=2.333335"
;#else
;  "cnst51=1.0"
;#endif

define loopcounter ncyc_max /* max value of ncyc used */
"ncyc_max=l8"
"DELTA8 = pwh_cp*2.0*cnst51*ncyc_max*2.0"

"l0=0"
"l1=0"
"l2=0"
"l3=0"
"l4=td1"
"l5=td2/8"


baseopt_echo
"acqt0=0"

aqseq 312



  ze 
1 d11 pl12:f2
2 100u do:f2

  20u pl11:f1
  (2mp ph1):f1
  20u
  (3mp ph2):f1

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  20u do:f1
  30u fq=0:f1

  d12 pl1:f1 pl2:f2
  50u UNBLKGRAD

  ; purge equilibrium magnetisation
  (p3 ph1):f2
  p16:gp2
  d16

2u
"ncyc_cp.idx=l1"
2u rpp11 rpp12 rpp13 rpp14 rpp21 rpp22 rpp23 rpp24

; l2 = ncyc
"l2 = (trunc(ncyc_cp + 0.3))"
2u
; l3 = ncyc_max - ncyc (unless ncyc = 0)
if "l2 > 0 " {
  "l3 = (trunc(ncyc_max - l2 + 0.3))"
}
else {
  "l3 = 0"
}

; calculate delay times
if "l2 > 0" {
  "tauCPMG = (time_T2*0.25)/l2"
  "tauCPMG1 = tauCPMG - (DELTA8*0.25 + 0.2u*l3)/l2"
}
else {
  "tauCPMG = time_T2*0.25"
  "tauCPMG1 = 2u"
}

  ; begin main sequence
  (p1 ph5):f1
  "DELTA1=d21-p1*0.6366"
  DELTA1

  (p3 ph6):f2
  if "l2 == 0" goto 51

  if "l0 % 4 == 0" goto 11
  if "l0 % 4 == 1" goto 21
  if "l0 % 4 == 2" goto 31
  if "l0 % 4 == 3" goto 41



11  (p4 ph1):f2 ; A
    d22
    ; begin CPMG (A)
; First half of CPMG
  if "l2 > 0" {
13   tauCPMG1
    ;(center (pwh_cp*2.0 ph11):f1 (p4 ph11):f2 )
    (center (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1 (p4 ph11):f2 )
    tauCPMG1 ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 13 times l2
  }
  if "l3 > 0" {
14   0.2u
    ;(center (pwh_cp*2.0 ph11):f1 (p4 ph11):f2 )
    (center (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1 (p4 ph11):f2 )
    0.2u ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 14 times l3

}
; Central 180o 1H and 13C pulses
  (center (pwh_cp*2.0 ph7):f1 (p4 ph2):f2 )
; Second half of CPMG
  if "l3 > 0" {
15   0.2u dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    ;(center (pwh_cp*2.0 ph21):f1 (p4 ph21):f2 )
    (center (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1 (p4 ph21):f2 )
    0.2u
    lo to 15 times l3
  }
  if "l2 > 0" {
16   tauCPMG1 dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    ;(center (pwh_cp*2.0 ph21):f1 (p4 ph21):f2 )
    (center (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1 (p4 ph21):f2 )
    tauCPMG1
    lo to 16 times l2
  }
    ;(center (p2 ph7):f1 (p4 ph1):f2)
    d22
    d0
    goto 99

21  d22         ; A'
; First half of CPMG
  if "l2 > 0" {
23   tauCPMG1
    (center (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1 (p4 ph11):f2 )
    ;(center (pwh_cp*2.0 ph11):f1 (p4 ph11):f2 )
    tauCPMG1 ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 23 times l2
  }
  if "l3 > 0" {
24   0.2u
    ;(center (pwh_cp*2.0 ph11):f1 (p4 ph11):f2 )
    (center (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1 (p4 ph11):f2 )
    0.2u ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 24 times l3

}
; Central 180o 1H and 13C pulses
  (center (pwh_cp*2.0 ph7):f1 (p4 ph2):f2 )
; Second half of CPMG
  if "l3 > 0" {
25   0.2u dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    ;(center (pwh_cp*2.0 ph21):f1 (p4 ph21):f2 )
    (center (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1 (p4 ph21):f2 )
    0.2u
    lo to 25 times l3
  }
  if "l2 > 0" {
26   tauCPMG1 dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    ;(center (pwh_cp*2.0 ph21):f1 (p4 ph21):f2 )
    (center (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1 (p4 ph21):f2 )
    tauCPMG1
    lo to 26 times l2
  }
    d22
    (p4 ph1):f2
    d0
    goto 99

31  d0           ; B
    d22
    ; begin CPMG (B)
; First half of CPMG
  if "l2 > 0" {
33   tauCPMG1
    ;(center (pwh_cp*2.0 ph11):f1 (p4 ph11):f2 )
    (center (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1 (p4 ph11):f2 )
    tauCPMG1 ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 33 times l2
  }
  if "l3 > 0" {
34   0.2u
    ;(center (pwh_cp*2.0 ph11):f1 (p4 ph11):f2 )
    (center (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1 (p4 ph11):f2 )
    0.2u ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 34 times l3

}
; Central 180o 1H and 13C pulses
  (center (pwh_cp*2.0 ph7):f1 (p4 ph2):f2 )
; Second half of CPMG
  if "l3 > 0" {
35   0.2u dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    ;(center (pwh_cp*2.0 ph21):f1 (p4 ph21):f2 )
    (center (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1 (p4 ph21):f2 )
    0.2u
    lo to 35 times l3
  }
  if "l2 > 0" {
36   tauCPMG1 dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    ;(center (pwh_cp*2.0 ph21):f1 (p4 ph21):f2 )
    (center (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1 (p4 ph21):f2 )
    tauCPMG1
    lo to 36 times l2
  }
    ;(center (p2 ph7):f1 (p4 ph1):f2)
    d22
    (p4 ph1):f2
    goto 99

41  d0          ; B'
    ; begin CPMG (B')
    (p4 ph1):f2
    d22
; First half of CPMG
  if "l2 > 0" {
43   tauCPMG1
    ;(center (pwh_cp*2.0 ph11):f1 (p4 ph11):f2 )
    (center (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1 (p4 ph11):f2 )
    tauCPMG1 ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 43 times l2
  }
  if "l3 > 0" {
44   0.2u
    ;(center (pwh_cp*2.0 ph11):f1 (p4 ph11):f2 )
    (center (pwh_cp ph12 pwh_cp*2.66667 ph11 pwh_cp ph12):f1 (p4 ph11):f2 )
    0.2u ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp23 ipp24
    lo to 44 times l3

}
; Central 180o 1H and 13C pulses
  (center (pwh_cp*2.0 ph7):f1 (p4 ph2):f2 )
; Second half of CPMG
  if "l3 > 0" {
45   0.2u dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    ;(center (pwh_cp*2.0 ph21):f1 (p4 ph21):f2 )
    (center (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1 (p4 ph21):f2 )
    0.2u
    lo to 45 times l3
  }
  if "l2 > 0" {
46   tauCPMG1 dpp11 dpp12 dpp13 dpp14 dpp21 dpp22 dpp23 dpp24
    ;(center (pwh_cp*2.0 ph21):f1 (p4 ph21):f2 )
    (center (pwh_cp ph22 pwh_cp*2.66667 ph21 pwh_cp ph22):f1 (p4 ph21):f2 )
    tauCPMG1
    lo to 46 times l2
  }
    d22
    goto 99



    ; ncyc = 0 - no relaxation delay
51  0.1u
    if "l0 % 4 == 0" goto 52
    if "l0 % 4 == 1" goto 53
    if "l0 % 4 == 2" goto 54
    if "l0 % 4 == 3" goto 55

52  (p4 ph1):f2 ; A
    d22
    (center (p2 ph7):f1 (p4 ph1):f2)
    d22
    d0
    goto 99

53  d22         ; A'
    (p4 ph1):f2
    (p2 ph7):f1
    d22
    (p4 ph1):f2
    d0
    goto 99

54  d0           ; B
    d22
    (center (p2 ph7):f1 (p4 ph1):f2)
    d22
    (p4 ph1):f2
    goto 99

55  d0          ; B'
    (p4 ph1):f2
    d22
    (p2 ph7):f1
    (p4 ph1):f2
    d22
    goto 99


99  (p3 ph1):f2
    "DELTA2=d21-4u"
    DELTA2 pl12:f2
    4u BLKGRAD

    go=2 ph31 cpd2:f2 
    20u do:f2

    d11 wr #0 if #0
    30u zd
    
    ; inner loop (pulse positions)
    30u iu0
    lo to 1 times 4

    ; middle loop (phases)
    30u ip6
    lo to 1 times 2

    ; loop over relaxation delays
    30u iu1
    lo to 1 times l4

    ; outer loop (d0)
    30u id0
    lo to 1 times l5

exit 

ph0=0
ph1=0
ph2=1
ph3=2
ph4=3
ph5=0 0 2 2
ph6=0 2
ph7=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph11=0 1 0 1
ph12=ph11-ph0
ph13=1 0 1 0
ph14=ph13-ph0
ph21=0 3 0 3
ph22=ph21+ph0
ph23=3 0 3 0
ph24=ph23+ph0
ph29=0
ph31=0 2 2 0 2 0 0 2


;pl3 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling (low power)
;p16: homospoil/gradient pulse                       [1 msec]
;p3: f2 channel - 90 degree high power pulse
;p1: f1 channel - 90 degree
;p2: f1 channel - 180 degree 
;d0 : incremented delay (2D) = in0/2-p3*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)CH
;d22 : 1/(8J)CH
;cnst4: = J(CH)
;inf1: 1/SW(C) = 2 * DW(C)
;in0: 1/ SW(C) = 2 * DW(C)
;nd0: 1
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEC
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62


;use gradient ratio:	gp 1 : gp 2
;			  11 :    7


;for z-only gradients:
;gpz1: 11%
;gpz2:  7%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100




;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1
