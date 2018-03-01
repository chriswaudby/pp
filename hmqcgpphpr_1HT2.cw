;methyl 1H T2 measurement
;option for (pseudo)1D measurement only (-DONE_D)
;L2 line (Tugarinov & Kay 2006)
;
;avance-version (07/04/04)
;HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive
;with decoupling during acquisition
;
;A. Bax, R.H. Griffey & B.L. Hawkins, J. Magn. Reson. 55, 301 (1983)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d3=1s/(cnst2*8)-p17-d16-larger(p1,p3)"
"d11=30m"
"d12=20u"
"d13=4u"

#ifndef ONE_D
"in0=inf2/2"
"d0=in0/2-0.63662*p3-2*p1"
#endif

"DELTA1=d2-p16-d16"
"DELTA2=d2-p16-d16-d12-4u-de+0.6366*p1"

define delay vdmin
"vdmin=4*(p1+p3+4u+p17+d16)"

"acqt0=de"
baseopt_echo

#ifndef ONE_D
aqseq 312
#endif

1 ze 
  vdmin
  d11 pl12:f2
2 d11 do:f2
  4u BLKGRAD

  20u
  "TAU1=vd*0.25-4u-p17-d16-p3"
  "TAU2=vd*0.25-p3-p1"
  "TAU3=vd*0.25-p1-4u-p17-d16-p3"
  "TAU4=vd*0.25-p3"

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  ; relaxation period
  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1 pl2:f2
  30u fq=0:f1
  50u UNBLKGRAD

  (p3 ph1):f2    ; crush eq'm 13C magnetisation
  d13
  p16:gp1
  d16

  ; start main sequence
  (p1 ph1):f1  ; INEPT
  DELTA1
  p16:gp2
  d16

  ; purge element
  (p3 ph11):f2
  d3
  p17:gp3
  d16
  (center (p2 ph1):f1 (p4 ph1):f2 )
  d3
  p17:gp3
  d16
  (p3 ph12):f2 

  ; t1 evolution
#ifndef ONE_D
  d0
#endif
  (p1 ph2):f1
  (p2 ph1):f1
  (p1 ph2):f1
#ifndef ONE_D
  d0
#endif
  (p3 ph13):f2

  ; relaxation period
  4u
  p17:gp4
  d16
  TAU1
  (p4 ph1):f2
  TAU2
  (p2 ph1):f1
  4u
  p17:gp4
  d16
  TAU3
  (p4 ph1):f2
  TAU4

  ; back-transfer
  d12 pl12:f2
  p16:gp2
  d16
  4u BLKGRAD
  DELTA2

  ; acquisition
  go=2 ph31 cpd2:f2 

#ifdef ONE_D
  d11 do:f2 mc #0 to 2
       F1QF(ivd)
#else
  d11 do:f2 mc #0 to 2
       F1QF(ivd)
       F2PH(ip13 & ip29, id0)
#endif /* ONE_D */
  4u BLKGRAD
exit 
  

ph1= 0 
ph2= 1 
ph11= 0 2
ph12= 1 1 3 3
ph13= 0
ph29=0
ph31=0 2


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse
;p17: gradient pulse [300 usec]
;p22 : f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d3 : 1/(8J)CH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz2: 7%
;gpz3: -40%
;gpz4: 29%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.32
;gpnam4: SMSQ10.32
                                          ;preprocessor-flags-start
;SINGLEDWELL: for initial sampling delay of one dwell-time with 
;	    option -DSINGLEDWELL (eda: ZGOPTNS)
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
;	    option -DOFFRES_PRESAT (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: hmqcphpr,v 1.4 2007/04/11 13:34:30 ber Exp $
