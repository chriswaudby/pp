;H(Z/D)QC
;for 13C, same processing as SFHZDQC
;run as pseudo-3D (td1 = 2), add and subtract to obtain Z/D components
;TODO check which is Z and D!
;
;Added option for off-resonance presat (e.g. to suppress urea signal), 21/6/15
;
;With option for 1D (first row)
;
;sfhmqcf2gpph
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

"p2=p1*2"

"in0=inf2"

"d0=in0/2-p3*4/3.1415"


"DELTA1=d21-1m-p16-d16-p1*0.6366"
"DELTA2=1m+p1*0.6366-de-4u"
"DELTA3=DELTA1-p2*0.5"
"acqt0=de"


"td1=2"
"l0=1"

aqseq 312



1 ze 
  d11 pl12:f2
2 10m do:f2

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  20u do:f1
  30u fq=0:f1

  d12 pl1:f1 pl2:f2
  50u UNBLKGRAD

  ; purge equilibrium magnetisation
  (p3 ph1):f2
  p16:gp2
  d16

  ; begin main sequence
  (p1 ph11):f1
  1m
  p16:gp1
  d16

  if "l0 %2 == 1"
     {
  (lalign (DELTA3 p2 ph13) (DELTA1 p3 ph12 d0 p3 ph1 DELTA1):f2 )
     }
  else
     {
  (ralign (p2 ph13 DELTA3 ) (DELTA1 p3 ph12 d0 p3 ph1 DELTA1):f2 )
     }

  DELTA2
  p16:gp1
  d16 pl12:f2
  4u BLKGRAD

  go=2 ph31 cpd2:f2 
  10m do:f2 mc #0 to 2 
     F1QF(ip12)
     F2EA(rp12 & iu0, id0)

exit 
  

ph1=0 
ph2=0 
ph11=0 
ph12=0 2 
ph13=0 0 1 1 2 2 3 3
ph29=0
ph31=0 2 2 0


;pl3 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling (low power)
;sp23: f1 channel - shaped pulse 120 degree 
;                   (Pc9_4_120.1000 or Q5.1000)
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p3: f2 channel -  90 degree high power pulse
;p1: f1 channel - 90 degree
;p2: f1 channel - 180 degree 
;d0 : incremented delay (2D) = in0/2-p3*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)CH
;cnst4: = J(CH)
;inf1: 1/SW(C) = 2 * DW(C)
;in0: 1/ SW(C) = 2 * DW(C)
;nd0: 1
;NS: 2 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEC
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62


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



