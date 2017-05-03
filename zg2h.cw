;zg2h
;avance-version (12/01/11)
;1D sequence
;using 2H lockswitch unit or BSMS 2H-TX board
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include<Avance.incl>
#include<Sysconf.incl>


"d11=30m"


"acqt0=-p1*2/3.1416"


1 ze
  d11 LOCKDEC_ON
  d11 LOCKH_ON
  d11 H2_PULSE

2 30m
  30m H2_LOCK 
  30m LOCKH_OFF
  d1
  30m LOCKH_ON
  d11 H2_PULSE 
  p1:D ph1
  go=2 ph31
  30m mc #0 to 2 F0(zd)

  d11 H2_LOCK
  d11 LOCKH_OFF
  d11 LOCKDEC_OFF
exit


ph1=0 2 2 0 1 3 3 1
ph2=0 2 2 0 1 3 3 1
ph3=2 0 0 2 3 1 1 3
ph31=0 2 2 0 1 3 3 1


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;ns: 1 * n, total number of scans: NS * TD0

;locnuc: off



;$Id: zg2h,v 1.14.8.1 2012/01/31 17:56:41 ber Exp $
