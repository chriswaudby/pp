; Dec 2016: use baseopt
;
; Jun 2013: added option for Ernst angle excitation
;
; Apr 2013: modified to use half-dwell first-point delay by default
;	    Added option for off-res presat
;
; Option for first row
;
; Removal of 13C equilibrium magnetisation (for methyl TROSY)
; Addition of clean-up gradient-pair
; Delays adjusted for zero first-order phase correction
; With options for 15N decoupling and 90,-180 or 180,-360 phase corr.
;
;hmqcphpr
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
"d2=1s/(cnst2*2)"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf1/2"

# ifndef ONE_D

#   ifdef LABEL_CN
    "p22=p21*2"
#   ifdef SINGLEDWELL
    "d0=in0-0.63662*p3-larger(p1,p21)"
#   else
    "d0=in0/2-0.63662*p3-larger(p1,p21)"
#   endif /*SINGLEDWELL*/
#   else
#   ifdef SINGLEDWELL
    "d0=in0-0.63662*p3-p1"
#   else
    "d0=in0/2-0.63662*p3-p1"
#   endif /*SINGLEDWELL*/
#   endif /*LABEL_CN*/

# endif /*ONE_D*/

# ifdef ERNST
    "p0=p1*cnst31"
# endif /*ERNST*/

"DELTA1=d2-p16-d16"
"DELTA2=d2-p16-d16-d12-4u-de+0.6366*p1"
"acqt0=de"

1 ze 
  d11 pl12:f2
2 d11 do:f2
  4u BLKGRAD

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1 pl2:f2
  30u fq=0:f1
  50u UNBLKGRAD
  (p3 ph1):f2
  d13
  p16:gp1
  d16*2 

# ifdef ERNST
  (p0 ph1):f1
# else
  (p1 ph1):f1
# endif /*ERNST*/

  DELTA1
  p16:gp2
  d16

# ifdef ONE_D

    ( center (p3 ph3 2u p3 ph4):f2 (p2 ph2):f1 )

# else

  (p3 ph3):f2
  d0

#   ifdef LABEL_CN
  (center (p2 ph5):f1 (p22 ph2):f3 )
#   else
  (p2 ph5):f1
#   endif /*LABEL_CN*/

  d0
  (p3 ph4):f2

# endif /*ONE_D*/

  d12 pl12:f2
  p16:gp2
  d16
  4u BLKGRAD
  DELTA2
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 F1PH(ip3 & ip29, id0)
  4u BLKGRAD
exit 
  
  
ph1=0 
ph2=0 
ph3=0 2
ph4=0 0 2 2
ph5=0 0 0 0 1 1 1 1
ph29=0

# ifdef FILTERED
ph31=0 0 0 0 2 2 2 2
# else
ph31=0 2 2 0 2 0 0 2
# endif /*FILTERED*/


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p16: homospoil/gradient pulse
;p22 : f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
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

;use gradient ratio:    gp 1 : gp 2
;                         31 :   7


;for z-only gradients:
;gpz1: 31%
;gpz2: 7%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
;SINGLEDWELL: for initial sampling delay of one dwell-time with 
;	    option -DSINGLEDWELL (eda: ZGOPTNS)
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
;	    option -DOFFRES_PRESAT (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: hmqcphpr,v 1.4 2007/04/11 13:34:30 ber Exp $
