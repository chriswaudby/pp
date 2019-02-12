;hsqcrexetf3gpsi3d.ts3.2.cw
;modifications for topspin 3.5pl7 (hopefully!)
;15N CW-CPMG HSQC
;recoded to use ncyc list instead of vdlist
;
;D. Flemming Hansen, Pramodh Vallurupalli, and Lewis E. Kay
;J. Phys. Chem. B, 2008, 112 (19), 5898-5904 
;
;----------------------------------------------------------
;CPMG frequencies specified via ncyc list (hard-coded in sequence below)
;va list, CW decoupling power for vCW = 2 * k * vCPMG 
;
;T_relax <= 40 ms (hard-coded in sequence below)
;
;warning - keep vCPMG <= 1000 Hz (2kHz 15N pulsing)
;
;----------------------------------------------------------

;prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


define delay TIME_T2
"TIME_T2 = 40m"

;vCPMG = ncyc / TIME_T2

;define list<loopcounter> ncyc = {0}
;define list<loopcounter> ncyc = {0 1 40}
define list<loopcounter> ncyc = {0 40 1 36 2 32 3 28 4 24 5 20 6 18 7 16 8 14 9 12 10}
define list<power> rflist=<$VALIST>

define delay DEL
"DEL=1m"



"p2=p1*2"
"p22=p21*2"
"p24=p25*0.5"
"d11=30m"
"d12=20u"
"d25=2.68m"
"d26=2.25m" 
"d28=5.0m"  ; tau_eq
"d17=p16"

"d10=6u"

"in10=inf2/2"


;"d31=0"  ; Feb 2014
"l11=0"  ; Feb 2014, loopcounter for cpmg


#   ifdef LABEL_CN
"DELTA1=d25-p16-d16+larger(p2,p8)+d10*2+10u"
#   else
"DELTA1=d25-p16-d16+p2+d10*2+10u"
#   endif /*LABEL_CN*/

"DELTA2=p16+d16+8u"
"DELTA4=p25*0.5-(p1*2/3.1416)" ;zeta
"DELTA5=(p45-4*p1)/3.1416" ;chi
"DELTA6=d26-p16-d16"
"DELTA7=d25-p16-d16"
"DELTA8=d28-p16-d16"
;"DELTA9=(1 / (1000*4) ) - p25/2000000" ;shortest CPMG delay, for temperature comp.

"spoff1=0"  ; H2O flip-down on-resonance
#ifdef LABEL_CN
"spoff13=bf2*((cnst21+cnst22)/2000000)-o2"  ; average of Ca/CO
#endif

"acqt0=0"    ; Feb 2014
baseopt_echo ; Feb 2014

aqseq 312


1 ze
  10m st0 ; interleaving
2 3m do:f3
3 3m do:f3
4 d11 do:f3

  4u pl16:f3
  4u fq=0:f1  ; 1H on H2O

  ; purge 1H magnetisation before d1
  4u pl12:f1
  (2mp ph1)
  (3mp ph2)
  4u pl1:f1

  ; recycle delay
  d1

  ; calculate delays etc. for CPMG elements
  4u
  "rflist.idx = l11"
if "ncyc==0" goto 71
  4u
  "DEL = TIME_T2/(4*ncyc) - p24"
71  4u


;  ; 15N compensation
;  ; initial delay to keep d1 constant:
;  if "l11==0" goto 74
;9   1000u  ; 2*(DELTA3+p25+DELTA3)
;    lo to 9 times COUNTER
;74  4u
;  ; now the 15N compensation pulses (unless d31 is at max 1000 Hz frequency)
;  if "d31 == 1000" goto 75
;     ; keep total number of pulses constant throughout experiment
;     ; (pulses applied at maximum frequency)
;     20u
;     "COUNTER2=d21*(1000-d31)*4 + 0.5"
;
;     4u pl23:f3
;     4u fq=cnst17(bf ppm):f3  ; 15N off-resonance
;
;8    DELTA9
;     (p25 ph2):f3
;     DELTA9 ;ipp20
;     lo to 8 times COUNTER2
;
;     4u fq=cnst16(bf ppm):f3  ; 15N on-resonance
;     4u ;rpp20
;75  10u

  ; crush Boltzmann Nz
  50u UNBLKGRAD
  4u pl1:f1 pl3:f3
  (p21 ph1):f3
  p16:gp0
  d16

  ; main sequence starts
  (p11:sp1 ph4:r):f1 ; flip-down
  2u pl1:f1

  (p1 ph1)
  p16:gp1
  d16
  DELTA6
  (center (p2 ph1) (p22 ph1):f3 )
  DELTA6
  p16:gp1
  d16
  (p1 ph2)

  p16:gp2
  d16 pl1:f1

  (p21 ph1):f3
  p16:gp3
  d16
  DELTA7
  (center (p2 ph1) (p22 ph1):f3 )
  DELTA7
  p16:gp3
  d16
  (p2 ph1) (p21 ph2):f3

; Nz equilibration
  p16:gp4
  d16 fq=cnst19(bf ppm):f1  ; 1H on amides
  DELTA8 pl23:f3  ; tau_eq (low power for 15N CPMG)

  (p1 ph5)
  DELTA5  ; chi
  (p1 ph1)
  DELTA4  ; zeta
  (p2 ph5)
  (p24 ph6):f3

  2u rflist:f1 

  ; first CPMG block
  if "ncyc==0" goto 81
     2u cpd1:f1 ph1
6    DEL
     (p25 ph2):f3
     DEL
     lo to 6 times ncyc
     2u do:f1
     goto 82
81   4u

  ; central refocusing pulse
82   2u pl33:f1 
     (center (p45 ph1) (p25 ph7):f3 )
     2u rflist:f1
  
  ; second CPMG block
  if "ncyc==0" goto 83
     2u cpd1:f1 ph1
7    DEL
     (p25 ph2):f3
     DEL
     lo to 7 times ncyc
     2u do:f1
     goto 84
83   4u
84   2u pl1:f1

  (p24 ph1):f3
  (p2 ph8)
  DELTA4 pl3:f3  ; zeta (back to high power on 15N)
  (p1 ph1)
  DELTA5    ; chi
  (p1 ph8)
  4u fq=0:f1  ; 1H on H2O

  ; re-equilibrate on Nz
  p16:gp5
  d16
  DELTA8  ; tau_eq

  ; back-transfer
  (p21 ph9):f3
  d10 gron6
  5u groff
  p16:gp7*EA
  d16
  DELTA7
# ifdef LABEL_CN
  (center (p2 ph1) (p8:sp13 ph1):f2 )
# else
  (p2 ph1)
# endif /*LABEL_CN*/
  d10 gron6*-1
  5u groff
  (p22 ph1):f3
  p16:gp7*-1*EA
  d16
  DELTA1   

  (center (p1 ph3) (p21 ph10):f3 )
  p16:gp8
  d16
  DELTA7
  (center (p2 ph1) (p22 ph1):f3 )
  DELTA7
  p16:gp8
  d16
  (center (p1 ph2) (p21 ph11):f3 )
  p16:gp9
  d16
  DELTA7
  (center (p2 ph1) (p22 ph1):f3 )
  DELTA7
  p16:gp9
  d16
  (p1 ph1)
  DELTA2
  (p2 ph1)
  4u
  p16:gp10
  d16 pl16:f3
  4u  BLKGRAD

;  go=2 ph31 cpd3:f3
;  d11 do:f3 mc #0 to 2
;    F1QF(calclc(l11, 1))
;    F2EA(calgrad(EA) & calph(ph11, +180), caldel(d10, in10) & calph(ph9, +180) & calph(ph31, +180))

  goscnp ph31 cpd3:f3
  3m do:f3

  if "ncyc > 0" goto 73
    ; reference plane - temperature compensation
     (p10:sp2 ph2:r):f1  ; H2O flip-down
     2u fq=cnst19(bf ppm):f1  ; 1H on amides
     2u pl33:f1 
     TIME_T2 cpd1:f1 ph1  ; 1H decoupling at pl33
     2u do:f1
     2u fq=0:f1  ; 1H on H2O
     (p10:sp3 ph4:r):f1  ; H2O flip-up
73   8u

  ; increment ncyc and loop counter for rflist
  3m st ncyc.inc iu11
  lo to 2 times nbl ; nbl = number of CPMG points

  ; reset loop counters and increment phase programs
  20u ncyc.res
  20u ru11
  3m ipp5 ipp6 ipp7 ipp8 ipp9 ipp31
  lo to 3 times ns  ; CPMG innermost loop, then ns

  d11 mc #0 to 4
     F1QF()
     F2EA(calgrad(EA) & calph(ph11, +180), caldel(d10, in10) & calph(ph9, +180) & calph(ph31, +180))



stop, exit


ph1=0 
ph2=1
ph3=2
ph4=3
ph5=1 1 3 3 
ph6=0 2 
ph7=0 0 2 2 
ph8=3 3 1 1
ph9=1 1 1 1 3 3 3 3 
ph10=0
ph11=1
ph31=0 2 0 2 2 0 2 0
  

;pl0 : 120dB

;1H pulses

;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;pl1 : f1 channel - power level for pulse (default)
;pl12 : f1 channel - power level for 10 kHz purge

;p10: f1 channel - 90 degree shaped pulse [1 ms]
;p11: f1 channel -  90 degree shaped pulse (7ms EBURP1 flipdown)
;sp1 : f1 channel - shaped pulse 90 degree (7ms EBURP1 flipdown)
;sp2 : f1 channel - 1 ms 90 degree flip-down
;sp3 : f1 channel - 1 ms 90 degree flip-back

;p45 : f1 channel - CW 180 degree pulse at pl33 (12-15 kHz)
;pl33 : f1 channel - CW decoupling power (12-15 kHz)
;valist : f1 channel - CW decoupling power (ca. 12-15 kHz as required for each ncyc point)


;15N pulses

;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;pl3 : f3 channel - power level for pulse (default)

;p24: f3 channel -  90 degree pulse at pl23 [45 usec]
;p25: f3 channel - 180 degree pulse at pl23 [90 usec]
;pl23: f3 channel - power level for spinlock

;pl16: f3 channel - power level for CPD/BB decoupling


;13C pulses

;p8 : f2 channel - 180 degree shaped pulse
;sp13: f2 channel - shaped pulse 180 degree  (Ca and C=O, adiabatic)


;gradients

;p16: homospoil/gradient pulse




;d1 : relaxation delay; 1-5 * T1
;d10 : incremented delay                             [6 usec]
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d17: gradient pulse length
;d21: length of CPMG mixing time (T_relax / 2)
;d24: 1/(4J)YH for YH
;     1/(8J)YH for all multiplicities
;d25: 1/(4J)YH for YH
;     1/(8J)YH for all multiplicities
;d26: 1/(4J(YH))
;d28: equilibration delay   [5 ms]
;cnst4: = J(YH)
;cnst16: o3p
;cnst17: o3p off-res. for compensation spinlock
;cnst19: centre of HN protons (offset, in ppm), e.g. 8.75 
;cnst11: for multiplicity selection = 4 for NH, 8 for all multiplicities
;cnst12: for multiplicity selection = 4 for NH, 8 for all multiplicities
;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)
;inf2: 1/SW(X) = 2 * DW(X)
;in10: 1/(2 * SW(X)) = DW(X)
;nd10: 2
;NS: 8 * n
;DS: >= 16
;td1: number of delays in vd-list
;td2: number of experiments in F2
;NBL: = td1
;FnMODE: QF in F1
;FnMODE: echo-antiecho in F2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence



;for z-only gradients:
;gpz0: 10%
;gpz1: 4%
;gpz2: -40%
;gpz3: -77.7%
;gpz4: 37.5%
;gpz5: -4.5%
;gpz6: 0.1%
;gpz7: 80%
;gpz8: 11.0%
;gpz9: 20.5%
;gpz10: 16.2%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam7: SMSQ10.100
;gpnam8: SMSQ10.100
;gpnam9: SMSQ10.100
;gpnam10: SMSQ10.100

;note: the values in the vd-list are interpreted as field-strength in Hz



;$Id: hsqcrexetf3gpsi3d,v 1.9 2009/07/17 16:37:47 ber Exp $
