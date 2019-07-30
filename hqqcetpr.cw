; HQQC (gradient-selected)
; Chris Waudby Jul 2019
;

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf1/2"
"d0=in0"

define list<grad_scalar> gl1 = { 0.13 }
define list<grad_scalar> gl2 = { 0.1 }
define list<grad_scalar> gl3 = { -0.87 0.73 }

"DELTA1=d2-p16-d16"
;"DELTA2=d2-p16-d16-d12-4u-de+0.6366*p1"
"acqt0=0.6366*p1"
baseopt_echo

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
  30u fq=0:f1

  d12 pl1:f1 pl2:f2
  50u UNBLKGRAD
  (p3 ph1):f2
  4u
  p16:gp1
  d16*2 

  (p1 ph1):f1

  ; note this is NOT an inept transfer - delays are 1/2J - transfer to CxHyHzHz
  DELTA1 ; 1/2J
  p16:gp2*gl1
  d16
  (center (p2 ph1):f1 (p3 ph1):f2 )
  p16:gp2*gl1
  d16
  DELTA1

  ; CxHyHzHz -> CxHyHxHx [4QC]
  p1 ph1
  d0
  p16:gp2*gl2
  d16
  (center (p2 ph1):f1 (p4 ph1):f2 )
  p16:gp2*gl2*-1
  d16
  p1 ph3

  p16:gp2*gl3
  d16
  DELTA1
  (center (p2 ph1):f1 (p3 ph1):f2 )
  p16:gp2*gl1
  d16 pl12:f2
  DELTA1 BLKGRAD

  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 F1EA(gl3.inc, id0)
  4u BLKGRAD
exit 
  
  
ph1=0 
ph3=2
ph29=0
ph31=0


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p16: homospoil/gradient pulse
;p22 : f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz2: 80%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100

