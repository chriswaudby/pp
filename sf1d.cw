;1D homonuclear SOFAST
;Added option for off-resonance presat (e.g. to suppress urea signal), 21/6/15
;
;With option for 1D (first row)
;
;sfhmqcf3gpph
;avance-version (09/11/18)
;SOFAST HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive
;with decoupling during acquisition
;
;P.Schanda and B. Brutscher, J. Am. Chem. Soc. 127, 8014 (2005)
;
;$CLASS=HighRes
;$DIM=2D
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
;"cnst39=0.529"

;"DELTA2=p39*cnst39-de-4u"
"DELTA2=p39*0.529-de-4u"
"acqt0=de"

# ifdef OFFRES_PRESAT
  "TAU=d1-10m-60u-d12*2-d13"
# else
  "TAU=d1-10m"
# endif /*OFFRES_PRESAT*/


;"spoff23=bf1*(cnst19/1000000)-o1"
;"spoff24=bf1*(cnst19/1000000)-o1"


1 ze 
  DELTA2
2 10m

# ifdef OFFRES_PRESAT
  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  ; TAU = d1
  TAU cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  30u fq=0:f1
# else
  ; TAU = d1
  TAU
# endif /*OFFRES_PRESAT*/

3 d12 
  50u UNBLKGRAD

  p16:gp2
  d16

  (p39:sp23 ph1):f1
  p16:gp1
  d16

  (p40:sp24 ph2):f1 

  DELTA2
  p16:gp1
  d16
  4u BLKGRAD

  go=2 ph31 
  10m mc #0 to 2 
exit 
  

ph1=0 2
ph2=0 0 1 1 2 2 3 3
ph29=0
ph31=0 2 2 0


;sp23: f1 channel - shaped pulse 120 degree (Pc9_4_120.1000)
;sp24: f1 channel - shaped pulse 180 degree (Reburp.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (2325 us at 800 MHz)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000              (1700 us at 800  MHz)
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;cnst19: 1H chemical shift for excitation (offset, in ppm)  [8.2 ppm]
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;NS: 2 * n
;DS: 16

;use gradient ratio:	gp 1 : gp 2
;			  11 :    7


;for z-only gradients:
;gpz1: 11%
;gpz2:  7%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


