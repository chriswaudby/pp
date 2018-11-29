;cpmghmqcph
;CPMG-HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum coherence
;xy16 CPMG blocks during coherence transfers
; turn on using -DCPMG1 and -DCPMG2
;phase sensitive
;with decoupling during acquisition

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p4=p3*2"
"d11=30m"
"d12=20u"
"d13=4u"

"d2=1s/(cnst2*2)"
"DELTA1=d2/32-larger(p1,p3)"
"DELTA=DELTA1*2"

"in0=inf1/2"
"d0=in0/2-p1-p3*0.6366"

"acqt0=0.6366*p1"
baseopt_echo



1 ze 
  4u pl1:f1 pl2:f2
2 d1 do:f2

  4u pl2:f2
  50u UNBLKGRAD
  (p3 ph1):f2    ; crush eq'm magnetisation
  20u
  p16:gp1
  d16
  4u BLKGRAD
  
  p1 ph1

#ifdef CPMG1
; phases: x y x y  y x y x  -x -y -x -y  -y -x -y -x
;         1 2 1 2  2 1 2 1  3  4  3  4   4  3  4  3
DELTA1
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA1
#else
  d2
#endif /* CPMG1 */

  p3:f2 ph13
  d0
  p2 ph12
  d0
  p3:f2 ph14

#ifdef CPMG2
; phases: x y x y  y x y x  -x -y -x -y  -y -x -y -x
;         1 2 1 2  2 1 2 1  3  4  3  4   4  3  4  3
DELTA1
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph2):f1 (p4 ph2):f2 ) ; y
DELTA
(center (p2 ph1):f1 (p4 ph1):f2 ) ; x
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA
(center (p2 ph4):f1 (p4 ph4):f2 ) ; -y
DELTA
(center (p2 ph3):f1 (p4 ph3):f2 ) ; -x
DELTA1 pl12:f2
#else
  d2 pl12:f2
#endif /* CPMG2 */

  go=2 ph31 cpd2:f2 
  ;d1 do:f2 mc #0 to 2 F1PH(calph(ph13, +90), caldel(d0, +in0))
  d1 do:f2 mc #0 to 2
      F1PH(calph(ph13, +90), caldel(d0, +in0) & calph(ph13, +180) & calph(ph31, +180))
exit 
  

ph1=0 
ph2=1
ph3=2
ph4=3
ph12=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph13=0 2
ph14=0 0 2 2 
ph31=0 2 2 0 2 0 0 2


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;d0 : incremented delay (2D)                  [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)XH
;d11: delay for disk I/O                             [30 msec]
;d13: short delay                                    [4 usec]
;cnst2: = J(XH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;ns: 16 * n
;ds: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence



;$Id: hmqcph,v 1.6 2012/01/31 17:49:23 ber Exp $
