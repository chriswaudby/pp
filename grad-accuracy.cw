;testing gradient +/- accuracy

#include <Avance.incl>
#include <Delay.incl>
#include <Grad.incl>


"d11=30m"


1 ze
2 d1 BLKGRAD
  50u UNBLKGRAD
  p1 ph1
  20u
  p16:gp1
  d16
  p16:gp2
  d16
  go=2 ph31
  wr #0
  d11 BLKGRAD
exit


ph1=0
ph31=0


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;p16: homospoil/gradient pulse [5 ms]
;gpz1: 50%
;gpz2: vary around -50%
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d16: delay for homospoil/gradient recovery

