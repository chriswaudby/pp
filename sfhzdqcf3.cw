;SOFAST-H(Z/D)QC
;run as pseudo-3D (td1 = 2), add and subtract to obtain Z/D components
;TODO check which is Z and D!
;
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
"d21=1s/(cnst4*2)"


"in0=inf2"

# ifdef ONE_D
"d0=6u"
#else
"d0=in0/2-p21*4/3.1415"
# endif /*ONE_D*/


"DELTA1=d21-p16-d16-p39*cnst39"
"DELTA2=p39*cnst39-de-4u"
"DELTA3=DELTA1-p40*0.5"
"acqt0=de"


# ifdef OFFRES_PRESAT
  "TAU=d1-10m-60u-d12*2-d13"
# else
  "TAU=d1-10m"
# endif /*OFFRES_PRESAT*/


"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"


"td1=2"
"l0=1"

aqseq 312



1 ze 
  d11 pl26:f3
2 10m do:f3

# ifdef OFFRES_PRESAT

  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  TAU cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  30u fq=0:f1

# else

  TAU

# endif /*OFFRES_PRESAT*/

3 d12 pl3:f3
  50u UNBLKGRAD

  ; purge Nz
  (p21 ph1):f3
  p16:gp2
  d16

  ; begin main sequence
  (p39:sp23 ph11):f1
  p16:gp1
  d16

  if "l0 %2 == 1"
     {
  (lalign (DELTA3 p40:sp24 ph1) (DELTA1 p21 ph12 d0 p21 ph1 DELTA1):f3 )
     }
  else
     {
  (ralign (p40:sp24 ph1 DELTA3 ) (DELTA1 p21 ph12 d0 p21 ph1 DELTA1):f3 )
     }

  DELTA2
  p16:gp1
  d16 pl26:f3
  4u BLKGRAD

  go=2 ph31 cpd3:f3 
  10m do:f3 mc #0 to 2 
     F1QF(ip12)
     F2EA(rp12 & iu0, id0)

exit 
  

ph1=0 
ph2=0 
ph11=0 
ph12=0 2 
ph31=0 2


;pl3 : f3 channel - power level for pulse (default)
;pl26: f3 channel - power level for CPD/BB decoupling (low power)
;sp23: f1 channel - shaped pulse 120 degree 
;                   (Pc9_4_120.1000 or Q5.1000)
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (3.0ms at 600.13 MHz)
;                  (or Q5.1000 (90o)            (2.0ms at 600.13 MHz) )
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000               (1.0ms at 600.13 MHz)
;d0 : incremented delay (2D) = in0/2-p21*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)NH
;cnst4: = J(NH)
;cnst19: H(N) chemical shift (offset, in ppm)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;           Q5.1000: -0.07
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/ SW(N) = 2 * DW(N)
;nd0: 1
;NS: 2 * n
;DS: 16
;aq: <= 50 msec
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEC
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence
;          use pulse of >= 350 usec


;use gradient ratio:	gp 1 : gp 2
;			  11 :    7


;for z-only gradients:
;gpz1: 11%
;gpz2:  7%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100




;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1



;$Id: sfhmqcf3gpph,v 1.1.2.8 2009/11/18 11:19:58 ber Exp $
