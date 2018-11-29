;3D methyl C(Z/D-quantum)CH SFH(Z/D)QC-NOESY-SFHMQC
;for methyl-methyl NOES
;run with F1 as QF, td = 4 * actual td (for interleaving)

;Option for NUS using Topspin 3
;Chris W, Oct 2018

;F1(H) -> F2(C[ZQ/DQ],t1) ---NOE--> F1(H) -> F2(C[mq],t2) -> F1(H,t3)
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



;------------options for first (in transfer pathway) 13C dim (F1)
"in0=inf1"		; first 13C dim 
"d0=in0/2-p3*4/3.1415"

;------------options for second (in transfer pathway) 13C dim (F2)
"in10=inf2"		; second 13C dim
"d10=in10/2-p3*4/3.1415"


; place pulses on-resonance
;"spoff23=bf1*(cnst19/1000000)-o1"
;"spoff24=bf1*(cnst19/1000000)-o1"
;"spoff25=bf1*(cnst19/1000000)-o1"

"TAU=d8-p16*2-d16*2-p3-8u"  ; noe mixing time


; delays for first SFHMQC
"DELTA1=d2-p41*cnst39-p16-d16"
"DELTA4=d2-p41*cnst39-p16-d16-p40*0.5"

; delays for second SFHMQC
"DELTA2=d2-p16-d16-p39*cnst39"
"DELTA3=p39*cnst39-de-4u"
"acqt0=de"

aqseq 312


1 ze
  d11 pl12:f2
  4u BLKGRAD
2 d11 do:f2
  d1 pl2:f2
  50u UNBLKGRAD

;-------------------------kill equm 13C magnetisation

  (p3 ph1):f2
  4u
  p16:gp1
  d16*2

;-------------------------start first 13C HMQC element
;uses DELTA1, DELTA4, DELTA5
  (p41:sp25 ph11):f1
  p16:gp2
  d16
  (lalign (DELTA4 p40:sp24 ph1) (DELTA1 p3 ph12 d0 p3 ph1 DELTA1):f2)
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


;------------------------start second 13C HMQC element (120o excitation)
;uses DELTA2, DELTA3
  (p39:sp23 ph1):f1
  p16:gp4
  d16

  (center (p40:sp24 ph2):f1 (DELTA2 p3 ph13 d10 p3 ph1 DELTA2):f2 )

  DELTA3
  p16:gp4
  d16 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2

  d11 do:f2 mc #0 to 2
    F1I(ip11, 2, ip12, 2)
    F1QF(id0)
    F2PH(rd0 & ip13, id10)
  4u BLKGRAD

exit


ph1= 0
ph2= 1
ph11=0
ph12=0 
ph13=0 2
ph31=0 2


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
;NS: 2 * n
;DS: 128
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI (or States) in F1
;FnMODE: States-TPPI (or States) in F2
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000  /  Reburp
;p41: f1 channel - 90 degree shaped pulse for excitation
;                      Pc9_4_90.1000 
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529


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
