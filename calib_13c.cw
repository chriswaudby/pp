; 13C 90 calibration
; set -DCALIB_C, adjust p20/pl20 for null, transmitter at cnst20 (ppm)
;
;hsqcphpr
;avance-version (07/04/04)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive
;with decoupling during acquisition
;
;G. Bodenhausen & D.J. Ruben, Chem. Phys. Lett. 69, 185 (1980)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Delay.incl>
#include <Grad.incl>


"p2=p1*2"
"d2=p2"
"p4=p3*2"
"p22=p21*2"
"d4=1s/(cnst2*4)"
"d11=30m"
"d12=20u"
"d13=4u"

"DELTA2=d4-larger(p2,p4)/2-p16-d16-4u"
"acqt0=p1*0.6366"

1 ze 
  d11 pl12:f2
2 d11 do:f2
3 d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1 pl2:f2 UNBLKGRAD
  p1 ph1
  4u
  p16:gp1
  d16
  DELTA2 
  (center (p2 ph2) (p4 ph5):f2 )
  DELTA2 
  p16:gp1
  d16
  4u
  (p1 ph3) 
  4u
  p16:gp2 ; zz purge
  d16

# ifdef CALIB_C
4u pl20:f2 fq=cnst20 (bf ppm):f2
(p20 ph1):f2
p16:gp2*1.3
d16
4u pl2:f2 fq=0:f2
# endif /* CALIB_C */

  (p3 ph6):f2
3u  
# ifdef LABEL_CN
  (center (p2 ph8) (p22 ph8):f3 )
# else
  (p2 ph8)
# endif /*LABEL_CN*/
3u
  (p3 ph7):f2
  4u
  p16:gp3
  d16
  (p1 ph4) 
  4u
  p16:gp4
  d16
  DELTA2 
  (center (p2 ph2) (p4 ph5):f2 )
  DELTA2 
  p16:gp4
  d16 BLKGRAD
  4u pl12:f2
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 F0(zd)
exit 
  

ph1=0
ph2=0
ph3=1
ph4=1
ph5=0 
ph6=0 2
ph7=0 0 0 0 2 2 2 2
ph8=0 0 2 2
ph29=0
ph31=0 2 0 2 2 0 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;sp3 : f2 channel - adiabatic inversion
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;cnst2: = J(XH)
;cnst20: 13C chemical shift of peak for calibration
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
;gpz1: 7 %
;gpz2: 50 %
;gpz3: 35 %
;gpz4: 13 %

;use gradient files:
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100

                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
;HALFDWELL: for initial sampling delay of half a dwell-time with
;           option -DHALFDWELL (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: hsqcphpr,v 1.4 2007/04/11 13:34:30 ber Exp $
