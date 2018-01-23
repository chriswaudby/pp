; pseudo-2D (indirect diffusion dimension, no carbon frequency dimension)
; With water saturation during diffusion delay
;
; Proton stimulated gradient-echo prior to HMQC
; H-C coupling evolution during encode/decode delays
; little delta limited to ~1.6 ms
;
; Removal of 13C equilibrium magnetisation (for methyl TROSY)
; Addition of clean-up gradient-pair
; Delays adjusted for zero first-order phase correction
; With options for 15N decoupling and 90,-180 or 180,-360 phase corr.
;
;hmqcphpr
;avance-version (07/04/04)
;HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive
;with decoupling during acquisition
;
;A. Bax, R.H. Griffey & B.L. Hawkins, J. Magn. Reson. 55, 301 (1983)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=

prosol relations=<triple_d>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<gradient> diff=<Difframp>

"p2=p1*2"
"d2=1s/(cnst2*2)"
"d11=30m"
"d12=20u"
"d13=4u"

"d0=3u"

"DELTA1=d2/2-p30-d16"
"DELTA2=d2-d12-4u-de+1.90986*p1"
"DELTA3=d20-d2/2-p19-d16-2*p1-d12-d13"

"acqt0=0"

1 ze 
2 d11 do:f2
  4u BLKGRAD
3 d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1 pl2:f2
  4u UNBLKGRAD
  (p3 ph1):f2
  d13
  p16:gp1
  d16 

  p1 ph6
  p30:gp6*diff
  d16
  DELTA1
  p1 ph1
  p19:gp2
  d16 pl9:f1
  DELTA3 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  p1 ph2
  p30:gp6*-1*diff
  d16
  DELTA1

  p3:f2 ph3
  d0

  (p2 ph5)

  d0
  p3:f2 ph4
  d12 pl12:f2
  4u BLKGRAD
  DELTA2
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 
	F1QF(igrad diff)
  4u BLKGRAD
exit 
  
  
ph1= 0 
ph2= 2 
ph3= 0 2
ph4= 0 0 0 0 2 2 2 2
ph5= 0 0 2 2 
ph6= 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
     1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3
ph29=0
ph31=0 2 0 2 2 0 2 0 2 0 2 0 0 2 0 2
     3 1 3 1 1 3 1 3 1 3 1 3 3 1 3 1


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p16: homospoil/gradient pulse
;p19: gradient pulse 2 (spoil gradient)
;p30: gradient pulse (little DELTA)
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;d20: diffusion time (big DELTA)
;cnst2: = J(CH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 8 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence

;use gradient ratio:    gp 1  : gp 2  : gp6
;                       -17.13: -13.17: var


;for z-only gradients:
;gpz1: -17.13%
;gpz2: -13.17%
;gpz6: 100%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam6: SMSQ10.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
;HALFDWELL: for initial sampling delay of half a dwell-time with 
;	    option -DHALFDWELL (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: hmqcphpr,v 1.4 2007/04/11 13:34:30 ber Exp $
