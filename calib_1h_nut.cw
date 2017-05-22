;calib_1h_nut.cw
;avance-version (02/05/31)
;1D sequence for 1H spinlock calibration by nutation of H2O
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Delay.incl>

define list<pulse> P1list = {0.1 15000 15500 16000 16100 16200 16300 16400 16500 16600 16700 16800 16900 17000 17100 17200 17500 18000 19000}

1 ze
2 30m
  d1 pl9:f1
  P1list ph1 
  go=2 ph31
  30m mc #0 to 2 
    F1QF(P1list.inc)
exit


ph1=0 2 2 0 1 3 3 1
ph31=0 2 2 0 1 3 3 1


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  high power pulse
;d1 : relaxation delay; 1-5 * T1
;NS: 1 * n, total number of scans: NS * TD0



;$Id: zg,v 1.8 2005/11/10 12:17:01 ber Exp $
