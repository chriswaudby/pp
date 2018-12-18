;cosyph
;avance-version (12/01/11)
;2D homonuclear shift correlation
;phase sensitive
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>


"in0=inf1"

"d0=3u"
"d2=p1*0.6366+3u"
"acqt0=-d2"

1 ze
2 d1
  p1 ph1
  d0
  go=2 ph31
  d1 mc #0 to 2 F1QF(id0) 
exit


ph1= 0 1 2 3
ph31=0 1 2 3


;pl1 : f1 channel - power level for pulse (default)
;p0 : f1 channel -  20 to 90 degree high power pulse
;p1 : f1 channel -  90 degree high power pulse
;d0 : incremented delay (2D)
;d1 : relaxation delay; 1-5 * T1
;inf1: 1/SW = 2 * DW
;in0: 1/(1 * SW) = 2 * DW
;nd0: 1
;ns: 4 * n
;ds: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ


;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1



;$Id: cosyph,v 1.9 2012/01/31 17:49:22 ber Exp $
