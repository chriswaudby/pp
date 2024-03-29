;sfhmqcf2gpph.nuws.cw
;1H,13C SOFAST HMQC
;switch on NUWS with -DNUWS
;option for first-increment only with -DONE_D 
;with apodisation-weighted sampling
;with exorcycle on 1H 180
;phase sensitive
;with decoupling during acquisition
;
;set sampling schedule via vclist:
;  total number of scans = ns * c
;  add one to first point of vclist to allow for dummy scans
;
;modified Chris Waudby 10/10/16
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

#ifdef NUWS
define loopcounter dsFlag
; dsFlag starts as 1, will be set to zero after first set of dummy scans is completed
"dsFlag=1"
#endif /* NUWS */

"d11=30m"
"d12=20u"
"d21=1s/(cnst2*2)"


"in0=inf1"

"d0=in0/2-p3*4/3.1415"


"DELTA1=d21-p16-d16-p39*cnst39"
"DELTA2=p39*cnst39-de-4u"
"acqt0=de"


"spoffs23=bf1*(cnst19/1000000)-o1"
"spoffs24=bf1*(cnst19/1000000)-o1"


1 ze
  d11 pl12:f2
2 d1 do:f2
3 d12 pl2:f2
  50u UNBLKGRAD

  ; purge equilibrium 13C
  (p3 ph1):f2
  4u
  p16:gp2
  d16

  (p39:sp23 ph1):f1
  p16:gp1
  d16

#ifdef ONE_D
  (center (p40:sp24 ph2):f1 (DELTA1 p3 ph3 0.1u p3 ph4 DELTA1):f2 )
#else
  (center (p40:sp24 ph2):f1 (DELTA1 p3 ph3 d0 p3 ph4 DELTA1):f2 )
#endif /* ONE_D */

  DELTA2
  p16:gp1
  d16 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2

#ifdef NUWS
  if "dsFlag==0" goto 10
  zd
  "dsFlag=0"
  goto 2
10 4u

  ; repeat acquisition block according to schedule in vclist
  lo to 2 times c
  30u ivc
#endif /* NUWS */

  d1 do:f2 mc #0 to 2
     F1PH(ip3, id0)
exit


ph1=0
ph2=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph3=0 2
ph4=0 0 2 2
ph31=0 2 2 0 2 0 0 2


;pl2 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling (low power)
;sp23: f1 channel - shaped pulse 120 degree (Pc9_4_120.1000)
;sp24: f1 channel - shaped pulse 180 degree (Reburp.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p3:  f2 channel -  90 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000
;d0 : incremented delay (2D) = in0/2-p3*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(CH) [125 Hz for methyls]
;cnst19: H(met) chemical shift (offset, in ppm)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;inf1: 1/SW(C) = 2 * DW(C)
;in0: 1/ SW(C) = 2 * DW(C)
;nd0: 1
;NS: 4 * n
;DS: 16
;aq: <= 50 msec
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEC
;cpd2: decoupling according to sequence defined by cpdprg2: garp4
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:	gp 1 : gp 2
;			  11 :    7


;for z-only gradients:
;gpz1: 11%
;gpz2:  7%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


                                          ;preprocessor-flags-start
                                          ;preprocessor-flags-end


;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1



;$Id: sfhmqcf3gpph,v 1.1.2.8 2009/11/18 11:19:58 ber Exp $
