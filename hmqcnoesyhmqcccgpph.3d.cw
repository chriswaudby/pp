;3D CCH HMQC-NOESY-13C HMQC
;for methyl-methyl NOES
;Option for NUS using Topspin 3
;Derived from hmqcnoesyhmqcccgpphpr.jk
;John K, Oct 2013
;Chris W, Dec 2016

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
"d11=45m"
"d12=20u"
"d13=4u"

"in0=inf1/2"		; first 13C dim
"in10=inf2/2"		; second 13C dim
;"in30=inf1/2"

;------------indirect 1H dim (F1)
;"d30=in30/2-p3"

;------------options for first (in transfer pathway) 13C dim (F2)
;"d0=in0/2-0.63662*p3-p1"
"d0=in0/2-0.63662*p3-p2"

;------------options for second (in transfer pathway) 13C dim (F3)
;"d10=in10/2-0.63662*p3-p1"
"d10=in10/2-0.63662*p3-p2"



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

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1 pl2:f2
  30u fq=0:f1
  50u UNBLKGRAD

;-------------------------kill equm 13C magnetisation

  (p3 ph1):f2
  d13
  p16:gp1
  d16*2

;-------------------------start first 13C HMQC element

  (p1 ph1):f1
  DELTA1
  p16:gp2
  d16

  ; 13C F2 evolution (MQ)
# ifdef NO_F1
;  ( center (p3 ph3 0.1u p3 ph4):f2 (p2 ph2):f1 )
  ( center (p3 ph3 0.1u p3 ph4):f2 (p1 ph1 p2 ph2 p1 ph1):f1 )
# else
  (p3 ph3):f2
  d0
;  (p2 ph2):f1
  (p1 ph1):f1
  (p2 ph2):f1
  (p1 ph1):f1
  d0
  (p3 ph4):f2
# endif /*NO_F1*/

  p16:gp2
  d16
  DELTA1

;------------------------start NOE period

  (p1 ph5):f1
  10u
  p16:gp3*0.71
  d16
  (p3 ph1):f2
  10u
  p16:gp3
  d16


  TAU


;------------------------start second 13C HMQC element

  (p1 ph6):f1

  DELTA3
  p16:gp4
  d16

# ifdef NO_F2
;  ( center (p3 ph7 0.1u p3 ph8):f2 (p2 ph2):f1 )
  ( center (p3 ph7 0.1u p3 ph8):f2 (p1 ph1 p2 ph2 p1 ph1):f1 )
# else
  (p3 ph7):f2
  d10
;  (p2 ph2):f1
  (p1 ph1):f1
  (p2 ph2):f1
  (p1 ph1):f1
  d10
  (p3 ph8):f2
# endif /*NO_F2*/

  d12 pl12:f2
  p16:gp4
  d16 
  4u BLKGRAD
  DELTA2
  go=2 ph31 cpd2:f2

  d11 do:f2 mc #0 to 2
	F2PH(ip7, id10)
	F1PH(rp7 & rd10 & ip3, id0)

  4u BLKGRAD

exit


ph1= 0
ph2= 1
ph3= 0 2
ph4= 0
ph5= 0
ph6= 0
ph7= 0 0 2 2
ph8= 0
ph29=0
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
;LABEL_CN: for 15N decoupling during indirect 13C evolution periods
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
;NOE_SAT: for water saturation during NOE mixing time
;F2_plane: for zero 13C phase evolution in F3
;F3_plane: for zero 13C phase evolution in F2
;F2_SINGLEDWELL: for single-dwell first-point delay in F2
;F3_SINGLEDWELL: for single-dwell first-point delay in F3
;TRIMPULSE: to apply trim-pulses in second 13C HMQC element, set p28
;NUS: for non-uniform sampling (Topspin 3)
                                          ;preprocessor-flags-end
