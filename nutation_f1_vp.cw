;nutation experiment on F1 channel using vplist
;for calibratoin of 1H B1 homogeneity
;avance-version (18/04/13)
;1D sequence
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=
;$RECOMMEND=y


#include <Avance.incl>
#include <Delay.incl>

define list<pulse> vplist = <$VPLIST>



1 ze
2 30m pl1:f1
  d1
  vplist
  go=2 ph31
  30m mc #0 to 2
    F1QF(vplist.inc)
exit


ph1=0
ph31=0


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;ns: 1 * n, total number of scans: NS * TD0



;$Id:$
