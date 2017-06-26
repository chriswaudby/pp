;zg with two selective flipdowns for solvent suppression
;avance-version (12/01/11)
;1D sequence
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=
;$RECOMMEND=y


#include <Avance.incl>
#include <Grad.incl>


"acqt0=-p1*2/3.1416"


1 ze
2 30m
  4u BLKGRAD
  d1
  20u UNBLKGRAD

(p11:sp11):f1
4u
p16:gp1
d16

(p12:sp12):f1
4u
p16:gp2
d16
4u BLKGRAD
4u pl1:f1

  p1 ph1
  go=2 ph31
  30m mc #0 to 2 F0(zd)
exit


ph1=0 2 2 0 1 3 3 1
ph31=0 2 2 0 1 3 3 1


;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;ns: 1 * n, total number of scans: NS * TD0



;$Id: zg30,v 1.12 2012/01/31 17:49:31 ber Exp $
