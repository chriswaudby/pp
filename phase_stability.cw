;testing phase stability
;avance-version (18/04/13)
;1D sequence
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=

#include <Avance.incl>
#include <Delay.incl>






1 ze
2 30m pl1:f1
  d1
  p1
  d20
  p1
  go=2 ph31
  30m mc #0 to 2
    F1QF()
exit


ph1=0
ph31=0


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;d20 : delay for phase evolution [1 ms]
;ns: 1, total number of scans: NS



;$Id:$
