;HMCM[CGCB]CACO.linear.2h.cw
; COSY experiment for assignment of methyl spin systems
;
; with selective 180 pulses
; explicit coding of t1 and t2 increments and quadrature
; frequency switching on 13C
;
;adapted by Chris Waudby Aug 2017
;
;3D sequence with
;   inverse correlation for triple resonance using inept transfer steps
;
;      F1(Hme) -> F2(Cme) -> F2(C->->Ca,t1)
;                          > F2(Ca->->C) -> F2(Cme,t2) -> F1(Hme,t3)
;
;phase sensitive (t1)
;phase sensitive (t2)
;using constant time in t2
;water suppression using watergate sequence
;
;V. Tugarinov & L.E. Kay, J. Am. Chem. Soc. 125, 13868-13878 (2003)
;K. Sinha, L. Jen-Jacobson & G.S. Rule. J. Biomol. NMR 56, 331�335 (2013)
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include<Avance.incl>
#include<Grad.incl>
#include<Delay.incl>




"p2=p1*2"
"p4=p3*2"

"d11=30m"  ; disk writing
"d12=20u"  ; short delay

"d4=1.8m"     ; 1/4J (methyl 1JCH = 125 Hz)
"DELTA1=d4-p16-d16-0.5*larger(p2,p4)"
"DELTA6=d4-p19-d16-0.5*larger(p2,p4)-8u-0.6366*p1-p11"
"DELTA7=d4-p19-d16-0.5*larger(p2,p4)-12u-p11"

"d24=7.0m"  ; 1/4J (1JCC = 35 Hz)
"DELTA2=d24-p14*0.5"
"d7=DELTA2-6u-p30"

"d23=5.75m"  ; optimal transfer time for CA->CO (with 52 Hz 1JCACO)
"DELTA3=d23-p14*1.5-1u"

"in0=inf2/2"  ; C' RT evolution
"d0=0.1u"

"in10=inf1/2"  ; Cmethyl CT evolution
"in30=in10"

"d10=in10/2-p3*0.6366-p1"  ; set initial delay for half-dwell
"d30=d24-in10/2-p3*0.6366-p14/2"
"DELTA5=d24-p14/2-p1-6u-p30"


"spoff1=0" ; 1H water selective pulses
"spoff2=0"
"spoff10=0"

"spoff3=0" ; 13C Q3 inversion, on-resonance for aliphatics
"spoff5=bf2*(cnst21/1000000)-o2" ; 13C' Q3 inversion, const21 = 173 ppm
"spoff7=o2-bf2*(cnst21/1000000)" ; 13CA Q3 inversion (C' on-res), const21 = 173 ppm
"spoff23=0" ; 13C' Q5 excitation, const21 = 173 ppm
"spoff25=0" ; 13C' Q5 excitation, const21 = 173 ppm

"p9=0.25*3.873/(bf2*(cnst21/1000000)-o2)"  ; 90 degree square pulse with null on C=O
;"p9=0.25*3.873/(bf2*(cnst21-cnst23)/1000000)"  ; 90 degree square pulse with null on C=O

"l1=td1/2"
"l2=td2/2"

aqseq 312

"acqt0=0"


1 ze

d11 LOCKDEC_ON
50u LOCKH_ON
d11 H2_PULSE
  d11 pl1:f1 pl12:f2 pl3:f3 pl17:f4
2 d11 do:f2
3 d12
4 d12
5 d12
  d11 H2_LOCK   ; lock-in over d1
  6m LOCKH_OFF

  ; relaxation delay
  d1

  ; turn lock hold on for 2H pulsing during main sequence
  d11 LOCKH_ON
  d11 H2_PULSE
  50u UNBLKGRAMP
  d12 pl1:f1 pl2:f2 pl3:f3 pl17:f4
  d12 ;fq=0:f2  ; 13C on methyls

  ;13C purge
  (p3 ph1):f2
  p16:gp1
  d16

  ;start sequence

  ;INEPT HM -> CM
  (p1 ph1)
  p16:gp2
  d16
  DELTA1
  (center (p2 ph1):f1 (p4 ph1):f2 )
  DELTA1
  p16:gp2
  d16 pl1:f1 pl2:f2
  (p1 ph2)

  (p11:sp1 ph3):f1  ; water flip-back
  p16:gp3
  d16 pl1:f1
  d12 ;fq=cnst23(bf ppm):f2  ; 13C on aliphatics

  ;INEPT CM -> CG
  (p3 ph11):f2
  DELTA2
  (p14:sp3 ph1):f2
  ;(p4 ph20):f2
  d7 pl2:f2
  ;switch on 2H decoupling
  (p30 ph2):f4
  3u
  3u cpds4:f4

#ifdef CB
  ;INEPT CG -> CB
  (p3 ph12):f2
  DELTA2
  (p14:sp3 ph1):f2
  DELTA2 pl2:f2
#ifdef CA
  ;INEPT CB -> CA
  (p3 ph12):f2
  DELTA2
  (p14:sp3 ph1):f2
  DELTA2 pl2:f2
#endif /* CA */
#endif /* CB */
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; BEGIN PROBLEMS!
  ; CA->CO
  (p3 ph12):f2
  1u
  (p14:sp5 ph1):f2  ; CO BS compensation
  DELTA3
  (p14:sp3 ph1):f2  ; CA
  1u
  (p14:sp5 ph1):f2  ; CO
  DELTA3 pl9:f2
  (p9 ph15):f2  ; Cali selective 90

  ;switch off 2H decoupling
  3u do:f4
  3u
  (p30 ph4):f4

  ; z-filter - do not move 13C carrier to C'
  3u fq=cnst21(bf ppm):f2
  p16:gp4
  d16

  ; t1 evolution (CO)
  (p13:sp23 ph16):f2  ; 90 CO selective
  d0 
  (center (p14:sp7 ph1):f2 (p22 ph1):f3 )  ; 180 on CA and 15N
  d0
;  (p15:sp6 ph1):f2  ; on-res 180
;  (p14:sp3 ph1):f2  ; on-res 180
;  3u
;  (p14:sp7 ph1):f2 ; 180 on CA (BS compensation)
;  3u 
  (p13:sp25 ph1):f2

  ; z-filter, move 13C carrier back to CA
  3u fq=0:f2
  p16:gp5
  d16

  ;switch on 2H decoupling
  (p30 ph2):f4
  3u pl9:f2
  3u cpd4:f4

  ; CO->CA 
  (p9 ph17):f2
  DELTA3
  (p14:sp5 ph1):f2  ; CO 
  1u
  (p14:sp3 ph1):f2  ; CA
  DELTA3 
  (p14:sp5 ph1):f2  ; CO BS compensation
  1u pl2:f2
;  DELTA4
;  (p14:sp5 ph1):f2  ; CO
;  3u
;  (p14:sp3 ph1):f2  ; CA
;  DELTA3 pl2:f2
  (p3 ph18):f2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; END PROBLEMS!

#ifdef CB
#ifdef CA
  ;INEPT CA -> CB [with refocussing of CG2]
  DELTA2
  (p14:sp3 ph1):f2
  DELTA2 pl2:f2
  (p3 ph18):f2
#endif /* CA */

  ;INEPT CB -> CG
  DELTA2
  (p14:sp3 ph1):f2
  DELTA2 pl2:f2
  (p3 ph18):f2
#endif /* CB */

  ;t2 evolution (with CT transfer CG -> CM)
  d10
  (p2 ph1):f1
  ;switch off 2H decoupling
  3u do:f4
  3u
  (p30 ph4):f4
  DELTA5
  (p14:sp3 ph1):f2
  d30 pl2:f2
  (p3 ph20):f2

  ; zz filter and water flip-down
  p16:gp6
  d16
  d12 ;fq=0:f2  ; 13C on methyls
  (p11:sp2 ph3):f1
  3u pl1:f1

  ;INEPT CM -> HM and watergate
  (p1 ph1)
  4u
  p19:gp7
  d16
  DELTA6
  (p11:sp10 ph3):f1
  4u pl1:f1
  (center (p2 ph1) (p4 ph1):f2 )
  4u
  (p11:sp10 ph3):f1
  4u
  DELTA7
  p19:gp7
  d16 pl12:f2
  4u BLKGRAMP

  go=2 ph31 cpd2:f2
  d11 do:f2 mc #0 to 2
    F1PH(ip20, id10 & dd30)
    F2PH(rp20 & rd10 & rd30 & ip16, id0)

  d11 H2_LOCK
  d11 LOCKH_OFF
  d11 LOCKDEC_OFF
exit


ph1 =0
ph2 =1
ph3 =2
ph4 =3
ph11=0 2
ph12=1 1 3 3
ph15=1
ph16=0 0 0 0 2 2 2 2
ph17=1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3
ph18=1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3
ph20=0
ph31=0 2 0 2 2 0 2 0 2 0 2 0 0 2 0 2


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;sp1: f1 channel - shaped pulse 90 degree (p11, H2O on resonance)
;sp2: f1 channel - shaped pulse 90 degree (p11, H2O on resonance)
;sp10: f1 channel - shaped pulse 90 degree (p11, H2O on resonance)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p11: f1 channel -  90 degree shaped pulse             [1 msec]


;pl2 : f2 channel - power level for pulse (default)
;pl9 : f2 channel - power level for p9
;pl12: f2 channel - power level for CPD/BB decoupling
;sp3: f2 channel - Q3.1000 CA refocusing
;sp5: f2 channel - Q3.1000 CO decoupling
;sp23: f2 channel - Q5.1000 90 degree excitation (C=O off resonance)
;sp25: f2 channel - Q5tr.1000 90 degree (C=O off resonance)
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p9 : f2 channel - 90 degree pulse with null on C=O
;p13: f2 channel - 90 degree shaped pulse
;p14: f2 channel - 180 degree shaped pulse [200 usec]
;p24: f2 channel - 180 degree shaped pulse [3000 usec]

;pl17: f4 channel - power level for 2H CPD decoupling
;p30: f4 channel - 90 degree pulse at pl17 (CPD 90)

;p16: homospoil/gradient pulse                         [1 msec]
;p19: gradient pulse 2                                 [500 usec]


;d0 : incremented delay (F1 in 3D)
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J(CH))                                       [2 msec]
;d10: incremented delay (F2 in 3D)
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d16: delay for homospoil/gradient recovery
;d24: 1/(4J(CmeC))                                     [7 msec]
;d25: 1/(4J'(CC))                                      [5 msec]
;d30: decremented delay (F2 in 3D)


;cnst23: Caliphatic chemical shift (30/40 ppm)
;o2p: Cmethyl chemical shift (19 ppm)
;in0: 1/(SW(Ca)) =  DW(Ca)
;nd0: 1
;in10: 1/(2 * SW(Cme)) = DW(Cme)
;nd10: 2
;in30: = in10
;NS: 8 * n [4 is possible]
;DS: >= 16
;td1: number of experiments in F1
;td2: number of experiments in F2       td2 max = 2 * d30 / in30
;FnMODE: States-TPPI (or TPPI) in F1
;FnMODE: States-TPPI (or TPPI) in F2
;cpd2: decoupling according to sequence defined by cpdprg2 [WALTZ-16]
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 67%
;gpz2: 11%
;gpz3: 50%
;gpz4: 40%
;gpz5: 30%


;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100



;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1

;PHC0(F2): 90
;PHC1(F2): -180
;FCOR(F2): 1

;SR(F1): (1/4) * SWH(F1)
;SR(F2): (1/4 + n) * SWH(F2), n = number of folding
