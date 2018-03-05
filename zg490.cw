;zg490.cw
;avance-version (02/05/31)
;1D sequence for pulse calibration
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


;$OWNER=ep
#include <Avance.incl>

"acqt0=-4*2/3.1416*p1"

1 ze
2 30m
  d1
;  p1 ph1 
;  3u
;  p1 ph1 
;  3u
;  p1 ph1 
;  3u
;  p1 ph1 
  p1*4 ph1
  go=2 ph31
  30m mc #0 to 2 F0(zd)
exit


ph1=0 2 2 0 1 3 3 1
ph31=0 2 2 0 1 3 3 1


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  high power pulse
;d1 : relaxation delay; 1-5 * T1
;NS: 1 * n, total number of scans: NS * TD0



;$Id: zg,v 1.8 2005/11/10 12:17:01 ber Exp $
