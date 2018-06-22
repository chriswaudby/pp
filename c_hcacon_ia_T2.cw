;c_hcacon_ia
; set td1 = 2*size of vdlist to allow for IPAP
; process with splitcomb (splitting in F1 = '0' option)
;avance-version (12/01/11)
;(HCa)CON
;2D sequence with
;   13C detected correlation for triple resonance using
;      multiple inept transfer steps
;
;      F2(H) -> F1(Ca) -> F1(C=O)
;            -> F3(N,t1) -> F1(C=O,t2)
;
;on/off resonance 13C pulses using shaped pulses
;using selective Ca pulse
;phase sensitive (t1)
;using IPAP scheme for virtual decoupling
;(use parameterset )
;
;W. Bermel, I. Bertini, I.C. Felli, R. Peruzzini & R. Pierattelli,
;   Chemphyschem. 11, 689-95 (2010)
;(W. Bermel, I. Bertini, V. Csizmok, I. C. Felli, R. Pierattelli & 
;   P. Tompa, J. Magn. Reson. 198, 275-281 (2009) )
;(;L. Duma, S. Hediger, A. Lesage & L. Emsley,
;   J. Magn. Reson. 164, 187-195 (2003) )
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple_c>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p4=p3*2"
"p22=p21*2"
"d11=30m"
"d12=20u"

"d3=1.1m"
"d4=1.8m"
"d22=4.5m"
"d23=12.4m"


"d0=3u"

"in0=inf2/2"


"DELTA=d0*2+p8"
"DELTA1=d22-d3-p4"
"DELTA2=d23-d22-p12"
"DELTA3=d23-p12-4u"
"DELTA4=d23/2-p12/2"


"l0=1"


"spoff13=bf1*((cnst22/2-cnst21/2)/1000000)"
"spoff23=0"
"spoff24=0"
"spoff25=0"
"spoff26=bf1*((cnst21-cnst22)/1000000)"
"spoff27=bf1*((cnst22-cnst21)/1000000)"
"spoff28=0"

"l6=td1/2"
"l3=td2/2"
aqseq 312


1 ze
  d11 pl12:f2 pl16:f3
2 d11 do:f2 do:f3
  d11
3 20u
4 20u
5 20u
6 20u
  d1 fq=cnst22(bf ppm):f1
  d12 pl2:f2 pl3:f3
  50u UNBLKGRAD

  (p3 ph1):f2
  d4 
  (center (p12:sp24 ph1) (p4 ph1):f2 )
  d4
  (p3 ph2):f2

  p16:gp1
  d16

  (p11:sp23 ph3)
  d3
  (p4 ph6):f2
  DELTA1
  (p12:sp26 ph6)
  4u
  (p25:sp28 ph1)
  d22
  (p12:sp26 ph1)
  4u
  (p11:sp25 ph1)

  p16:gp2
  d16 fq=cnst21(bf ppm):f1

  (p11:sp23 ph5)
  d22
  (p12:sp27 ph1)
  DELTA2
  (center (p12:sp24 ph1) (p22 ph1):f3 )
  DELTA3
  (p12:sp27 ph1)
  4u
  (p11:sp25 ph1)

  p16:gp3
  d16

  (p21 ph4):f3
  d0
  (center (p8:sp13 ph6) (p4 ph6):f2 )
  d0
  (p22 ph1):f3
  DELTA
  (p21 ph1):f3

  p16:gp4
  d16 

  if "l0 %2 == 1"
     {
     (p11:sp23 ph1)
     DELTA4
     (p12:sp27 ph1)
     DELTA4
     (center (p12:sp24 ph1) (p22 ph1):f3 )
     DELTA4
     (p12:sp27 ph1)
     DELTA4 ;pl16:f3
     }
  else
     {
     (p11:sp23 ph7)
     d22
     (p12:sp27 ph1)
     DELTA2
     (center (p12:sp24 ph1) (p22 ph1):f3 )
     DELTA4
     DELTA4 ;pl16:f3
     (p12:sp27 ph1)
     }

  ; relaxation period
  vd*0.25
  (center (p12:sp27 ph1) (p22 ph1):f3 (p4 ph1):f2)
  vd*0.25
  (p12:sp24 ph8)
  vd*0.25
  (center (p12:sp27 ph1) (p22 ph1):f3 (p4 ph1):f2)
  vd*0.25 pl12:f2 pl16:f3

  4u BLKGRAD

  go=2 ph31 cpd2:f2 cpd3:f3
  d11 do:f2 do:f3
  d11 wr #0 if #0 zd

  ; inner loop - IPAP
  20u iu0
  lo to 3 times 2

  ; middle loop - vdlist
  20u ivd
  lo to 4 times l6

  ; outer loop - states-tppi
  20u ip4
  lo to 5 times 2

  20u id0
  lo to 6 times l3


;  go=2 ph31 cpd2:f2 cpd3:f3
;  d11 do:f2 do:f3 mc #0 to 2 
;     F1QF(ivd)
;     F2I(iu0, 2)
;     F2PH(calph(ph4, +90), caldel(d0, +in0)) 
exit


ph1=0
ph2=1
ph3=0 0 0 0 2 2 2 2
ph4=0 0 2 2
ph5=0 2            
ph6=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph7=3
ph8=0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
ph31=0 2 2 0 2 0 0 2 0 2 2 0 2 0 0 2
     2 0 0 2 0 2 2 0 2 0 0 2 0 2 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;pl16: f3 channel - power level for CPD/BB decoupling
;sp13: f1 channel - shaped pulse 180 degree  (adiabatic)
;sp23: f1 channel - shaped pulse  90 degree  (on resonance)
;sp24: f1 channel - shaped pulse 180 degree  (on resonance)
;sp25: f1 channel - shaped pulse  90 degree  (on resonance)
;                   for time reversed pulse
;sp26: f1 channel - shaped pulse 180 degree  (C=O off resonance)
;sp27: f1 channel - shaped pulse 180 degree  (Ca off resonance)
;sp28: f1 channel - shaped pulse 180 degree  (Ca on resonance)
;                   sp28 requires higher selectivity than sp24
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p8 : f1 channel - 180 degree shaped pulse for inversion (adiabatic)
;p11: f1 channel -  90 degree shaped pulse
;p12: f1 channel - 180 degree shaped pulse
;p16: homospoil/gradient pulse                         [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p25: f1 channel - 180 degree shaped pulse (Ca, sp28)
;d0 : incremented delay (F1 in 2D)                     [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d3 : 1/(6J(HCa))                                      [1.1 msec]
;d4 : 1/(4J(HCa))                                      [1.8 msec]
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d16: delay for homospoil/gradient recovery
;d22: 1/(4J(COCa))                                     [4.5 msec]
;d23: 1/(4J(NCO))                                      [12.4 msec]
;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)
;o1p: CO chemical shift (cnst21)
;l0: flag to switch between inphase and antiphase
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(N)) = DW(N)
;nd0: 2
;ns: 16 * n
;ds: >= 32
;td1: number of experiments in F1 * 2
;FnMODE: States-TPPI (or TPPI) in F1
;cpd2: decoupling according to sequence defined by cpdprg2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 50%
;gpz2: 30%
;gpz3: 19%
;gpz4: 13%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100



;use AU-program splitcomb [ipap 2] to process data



;$Id: c_hcacon_ia,v 1.2 2012/01/31 17:49:22 ber Exp $
