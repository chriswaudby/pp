;4D HCCH SFHMQC-NOESY-SFHMQC
;for methyl-methyl NOES
;Option for NUS using Topspin 3
;Derived from hmqcnoesyhmqcccgpphpr.jk
;Chris W, July 2018

;F1(H) -> F2(C[mq],t1) ---NOE--> F1(H) -> F2(C[mq],t2) -> F1(H,t3)
;
;Indirect evolution order is t2, t1 (13Cdir, 13Cnoe)
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
"d11=45m" ; for disk access


;------------indirect 1H dim (F1)
"in0=inf1/2"
"d0=in0/2-p3"

;------------options for first (in transfer pathway) 13C dim (F2)
"in10=inf2/2"		; first 13C dim (NB two d10 delays present)
"d10=in10/2-p3*4/3.1415"

;------------options for third (in transfer pathway) 13C dim (F3)
"in30=inf3"		; second 13C dim (NB only one d30 delay present)
"d30=in30/2-p3*4/3.1415"


"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"
"spoff25=bf1*(cnst19/1000000)-o1"

"TAU=d8-p16*2-d16*2-p3-8u"  ; noe mixing time


; delays for first SFHMQC
"d5=d2-p39*cnst39-4u-p16-d16"
"DELTA1=d5+p3+p4-p40*0.5"

; delays for second SFHMQC
"DELTA2=d2-p16-d16-p39*cnst39"
"DELTA3=p39*cnst39-de-4u"
"acqt0=de"

;aqseq 4321	; for info only


1 ze
  d11 pl12:f2
  4u BLKGRAD
2 d11 do:f2
  d1
  50u UNBLKGRAD

;-------------------------kill equm 13C magnetisation

  (p3 ph1):f2
  4u
  p16:gp1
  d16*2

;-------------------------start first 13C HMQC element

  (p41:sp25 ph11):f1
  4u
  p16:gp2
  d16

  ; 1H F1 and 13C F2 evolution (MQ)
  (lalign
    (DELTA1 2*d0 d10 p40:sp24 ph2):f1
    (d5 p3 ph12 d0 p4 ph1 d0 2*d10 p3 ph1 d5):f2
  )

  4u
  p16:gp2
  d16
  (p41:sp25 ph1):f1

;------------------------start NOE period
  4u
  p16:gp3*0.71
  d16
  (p3 ph1):f2
  4u
  p16:gp3
  d16


  TAU


;------------------------start second 13C HMQC element
  (p39:sp23 ph1):f1
  p16:gp4
  d16

  (center (p40:sp24 ph2):f1 (DELTA2 p3 ph13 d30 p3 ph1 DELTA2):f2 )

  DELTA3
  p16:gp4
  d16 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2

#ifndef NUWS
  d11 do:f2 mc #0 to 2
    F3PH(ip13, id30)
	F2PH(rp13 & rd30 & ip12, id10)
	F1PH(rp13 & rd30 & rp12 & rd10 & ip11, id0)
#else
  d11 do:f2 mc #0 to 2
	F1PH(calph(ph11), caldel(d0))
	F2PH(calph(ph12), caldel(d10))
	F3PH(calph(ph13), caldel(d30))
#endif /*NUS*/
  4u BLKGRAD

exit


ph1= 0
ph2= 1
ph11=0
ph12=0 2
ph13=0 0 2 2
ph31=0 2 2 0


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
;p22: f3 channel - 180 degree high power pulse
;p28: f1 channel - trim pulse                        [1 msec]
;d0 : incremented delay (4D)
;d10: incremented delay (4D)
;d20: decremented delay (4D)
;d28: incremented delay (4D)
;d30: incremented delay (4D)
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d8 : mixing time
;d11: delay for disk I/O                             [45 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;inf1: 1/SW(H)
;inf2: 1/SW(C) = 2 * DW(C)
;inf3: 1/SW(C) = 2 * DW(C)
;in0: 1/(2 * SW(C)) = DW(C)
;in10: 1/(2 * SW(C)) = DW(C)
;in30: 1/(2 * SW(H)) = DW(H)
;nd0: 2
;nd10: 2
;NS: 4 * n
;DS: 16
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI (or States) in F1
;FnMODE: States-TPPI (or States) in F2
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz2: 29%
;gpz3: 23%
;gpz4: 13%
;gpz5: 43%
;gpz6: 19%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SINE.50
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam5: SINE.50
;gpnam6: SINE.100

                                          ;preprocessor-flags-start
;NUS: for non-uniform sampling (Topspin 3)
                                          ;preprocessor-flags-end
