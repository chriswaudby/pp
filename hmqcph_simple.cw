;hmqcph
;avance-version (12/01/11)
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


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"d11=30m"
"d12=20u"
"d13=4u"

"d2=1s/(cnst2*2)"
"DELTA1=d2-p16-d16"
"DELTA2=d2-p16-d16-d12-4u-de+0.6366*p1"

"in0=inf1/2"
"d0=in0/2-p1-p3*0.6366"

"acqt0=de"
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
  
  p1 ph1

  DELTA1
  p16:gp2
  d16

  p3:f2 ph3
  d0
  p2 ph2
  d0
  p3:f2 ph4

  d12 pl12:f2
  p16:gp2
  d16
  4u BLKGRAD
  DELTA2

  go=2 ph31 cpd2:f2 
  d1 do:f2 mc #0 to 2 F1PH(calph(ph3, +90), caldel(d0, +in0))
exit 
  

ph1=0 
ph2=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph3=0 2
ph4=0 0 2 2 
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
;ns: 4 * n
;ds: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence



;$Id: hmqcph,v 1.6 2012/01/31 17:49:23 ber Exp $
