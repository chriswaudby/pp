;1D SOFAST with pre-sat for acquisition of both 15N-edited and 
;15N-filtered 1Ds
;***** Use 90deg PC-9 excitation pulse *****

;sfhmqcf3gpph
;avance-version (08/03/03)
;SOFAST HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive
;with decoupling during acquisition
;
;P.Schanda and B. Brutscher, J. Am. Chem. Soc. 127, 8014 (2005)
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"d11=30m"
"d12=20u"
"d13=4u"
"d21=1s/(cnst4*2)"


"DELTA1=d21-p16-d16-p41*cnst41"
"DELTA2=d21-p16-d16-de-4u"

"spoff24=bf1*(cnst19/1000000)-o1"
"spoff25=bf1*(cnst19/1000000)-o1"


1 ze 
  d11 pl26:f3
2 d11 do:f3
3 d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl0:f1 pl3:f3
  50u UNBLKGRAD
  
  p16:gp2
  d16

  (p41:sp25 ph1):f1
  p16:gp1
  d16
  DELTA1

  (center (p40:sp24 ph2):f1 (p21 ph3 4u p21 ph4):f3 )

  DELTA2 
  p16:gp1
  d16 pl26:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3 
  d11 do:f3 mc #0 to 2 F0(zd)
exit 
  

ph1=0 0 0 0 2 2 2 2 1 1 1 1 3 3 3 3
ph2=0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
ph3=0 2
ph4=0 0 2 2 
ph29=0
# ifdef FILTERED
ph31=0 0 0 0 2 2 2 2 3 3 3 3 1 1 1 1
     0 0 0 0 2 2 2 2 3 3 3 3 1 1 1 1
     2 2 2 2 0 0 0 0 1 1 1 1 3 3 3 3
     2 2 2 2 0 0 0 0 1 1 1 1 3 3 3 3
# else
ph31=0 2 2 0 2 0 0 2 3 1 1 3 1 3 3 1
     0 2 2 0 2 0 0 2 3 1 1 3 1 3 3 1
     2 0 0 2 0 2 2 0 1 3 3 1 3 1 1 3
     2 0 0 2 0 2 2 0 1 3 3 1 3 1 1 3
# endif /*FILTERED*/


;pl3 : f3 channel - power level for pulse (default)
;pl26: f3 channel - power level for CPD/BB decoupling (low power)
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;sp25: f1 channel - shaped pulse 90 degree 
;                   (Pc9_4_90.1000)
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse                       [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p41: f1 channel -  90 degree shaped pulse for excitation
;                      Pc9_4_90.1000 (90o)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)NH
;cnst4: = J(NH)
;cnst19: H(N) chemical shift (offset, in ppm)
;cnst41: compensation of chemical shift evolution during p41
;           Pc9_4_90.1000: 0.514
;NS: 2 * n
;DS: 16
;aq: <= 50 msec
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence
;          use pulse of >= 350 usec


;use gradient ratio:	gp 1 : gp 2
;			  11 :   37


;for z-only gradients:
;gpz1: 11%
;gpz2: 37%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
