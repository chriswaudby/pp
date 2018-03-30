;3D HCH SFHMQC-NOESY-13C HMQC
;for methyl-methyl NOES
;for perdeuterated proteins in D2O
;selective excitation in first block; hard pulses in second
;Chris W, Mar 2018
;
;adapted from Rossi JBNMR 2016 Fig 4a
;
;at 800 MHz, with 0.7ppm 1H offset:
; p39 (Pc9.90): 2006 us
; p40 (Reburp): 1100 us

;F1(H,t1) ---NOE--> F1(H) -> F2(C[mq],t2) -> F1(H,t3)
;
;Indirect evolution order is t2, t1 (13C, 1Hnoe)
;Uses half-dwell first-point delay by default in all indirect dims
;Option for off-res presat
;Options for 2D planes in each 13C dim
;  (set both to get 1D or 2D HH plane with no 13C phase evolution)
;Removal of 13C equilibrium magnetisation
;Delays adjusted for zero first-order phase correction in acqusition dim


;$CLASS=HighRes
;$DIM=4D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d11=45m"
"d12=20u"
"d13=4u"


;------------options for first (in transfer pathway) 1H noe (F1)
"in0=inf1"		; first 1H dim
"d0=3u"

;------------options for second (in transfer pathway) 13C (F2)
"in10=inf2"		; 13C dim
"d10=in10/2-1.27324*p3"

;------------shaped pulses
"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"
"spoff25=bf1*(cnst19/1000000)-o1"
"spoal23=1"
"spoal24=0.5"
"spoal25=0"

;------------delays
"TAU=d8-p16*2-d16*2-p3-20u"
"DELTA1=d2-p16-d16"
"DELTA2=d2-p16-d16-d12-4u-de+0.63662*p1"
"DELTA3=d2-p16-d16"
"acqt0=de"

aqseq 321	; for info only


1 ze
  d11 pl12:f2
2 d11 do:f2
  4u BLKGRAD

#ifdef OFFRES_PRESAT
  30u fq=cnst21(bf hz):f1
#endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1 
  30u fq=0:f1
  50u UNBLKGRAD

;-------------------------initial purge gradient
  p16:gp1
  d16

;-------------------------1H evolution (t1)
  d16 pl12:f2
  0.1u cpd2:f2  ; 13C decoupling on

  (p39:sp23 ph11):f1 
  d0
  (p40:sp24 ph13):f1
  3u
  (p39:sp25 ph1):f1
  0.1u do:f2
  4u pl2:f2 
;------------------------start NOE period
  p16:gp3*0.71
  d16
  (p3 ph1):f2
  10u
  p16:gp3
  d16
  TAU
;------------------------start second 13C HMQC element
  4u pl1:f1 
  (p1 ph1):f1
  DELTA3
  p16:gp4
  d16
  ( center (p3 ph12 d10 p3 ph1):f2 (p2 ph1):f1 )
  d12 pl12:f2
  p16:gp4
  d16 
  4u BLKGRAD
  DELTA2
  go=2 ph31 cpd2:f2
  d11 do:f2 mc #0 to 2
	F2PH(ip12, id10)
	F1PH(rp12 & rd10 & ip11, id0)

  4u BLKGRAD

exit


ph1= 0
ph2= 1
ph11=0 2
ph12=0 0 2 2
ph13=0 0 0 0 1 1 1 1
ph29=0
ph31=0 2 2 0 2 0 0 2


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse                       [1 msec]
;p19: second gradient pulse                          [250 usec]
;d0 : incremented delay (1H)
;d10: incremented delay (13C)
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d8 : mixing time
;d11: delay for disk I/O                             [45 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;sp23: f1 channel - shaped pulse 90 degree (Pc9_4_90.1000)
;spoal23: 1
;sp25: f1 channel - shaped pulse 90 degree (Pc9_4_90.1000)
;spoal25: 0
;sp24: f1 channel - shaped pulse 180 degree (Reburp.1000)
;cnst19: H(methyl) offset for selective excitation, in ppm
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;NS: 8 * n
;DS: 16
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI (or States) in F1
;FnMODE: States-TPPI (or States) in F2
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: -31%
;gpz3: 43%
;gpz4: -23%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100

                                          ;preprocessor-flags-start
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
                                          ;preprocessor-flags-end
