; With option for off-res presat (removed orthogonal phase pre-sat)
;
;zgpr
;avance-version (06/11/09)
;1D sequence with f1 presaturation
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Delay.incl>


"d12=20u"


"acqt0=-p1*2/3.1416"


1 ze
2 10m

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  4u do:f1
  d12 pl1:f1
  30u fq=0:f1
  p1 ph1
  go=2 ph31
  10m mc #0 to 2 F0(zd)
exit


ph1=0 2 2 0 1 3 3 1
ph29=0
ph31=0 2 2 0 1 3 3 1


;pl1 : f1 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;p1 : f1 channel -  90 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;d12: delay for power switching                      [20 usec]
;NS: 1 * n, total number of scans: NS * TD0



;$Id: zgpr,v 1.9 2006/11/10 10:56:44 ber Exp $
