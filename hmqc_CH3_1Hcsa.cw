;IPAP HMQC for measurement of 1H CSA via 1H CSA/1H-13C DD CCR
; set td1 = 2*number of relaxation time points
;Chris Waudby, July 2020
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
"d3=1s/(cnst2*8)"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf2/2"
"d0=in0/2-0.63662*p3-2*p1"


; loop counter for IPAP
"l1=0"

define delay vdmin
"vdmin=2*(p2+4u+p17+d16)"

"acqt0=0"
baseopt_echo

aqseq 312

1 ze 
  vdmin
  d11 ph1:f1 pl3:f2
2 d11 

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
  "DELTA1=d2-p16-d16+0.6366*p1"
  DELTA1
  p16:gp2
  d16

  ; purge element
  (p3 ph11):f2
  "DELTA=d3-p17-d16-larger(p1,p3)"
  DELTA
  p17:gp3
  d16
  (center (p2 ph1):f1 (p4 ph1):f2 )
  DELTA
  p17:gp3
  d16
  (p3 ph12):f2 

  ; t1 evolution
  d0
  (p1 ph13):f1
  (p2 ph14):f1
  (p1 ph13):f1
  d0
  (p3 ph15):f2

  ; relaxation period
  "TAU=vd*0.5-p17-d16-p2-4u"
  TAU
  4u
  p17:gp4
  d16
  (p1 ph1):f1
  (p2 ph2):f1
  (p1 ph1):f1
  4u
  p17:gp4
  d16
  TAU
 

  ; IPAP back-transfer
  if "l1 % 2 == 0" {
  ; IP
  "DELTA2=d2*0.5-p16-d16-p3"
  p16:gp2
  d16
  DELTA2
  p4:f2 ph1
  "DELTA3=d2*0.5-p3-4u"
  DELTA3
  4u BLKGRAD
  } else {
  ; AP
  "DELTA2=d2-p16-d16-p3-4u"
  p16:gp2
  d16
  DELTA2
  p3:f2 ph1
  4u BLKGRAD
  }

  ; acquisition
  go=2 ph31
  d11 mc #0 to 2
    F1I(iu1, 2)
    F1QF(ivd)
    F2PH(ip15, id0)

  4u BLKGRAD
exit 
  

ph1= 0 
ph2= 1 
ph11=0 2
ph12=1 1 1 1 3 3 3 3
ph13=0 0 1 1 2 2 3 3
ph14=1 1 2 2 3 3 0 0
ph15=0
ph29=0
ph31=0 2 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
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
;NS: 8 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ


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
