; Nov 2014: L1/L3 purge element (Korzhnev 2004)
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
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d3=1s/(cnst2*8)-p17-d16-larger(p1,p3)"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf1/2"

#   ifdef SINGLEDWELL
    "d0=in0-0.63662*p3-p1"
#   else
    "d0=in0/2-0.63662*p3-p1"
#   endif /*SINGLEDWELL*/


"DELTA1=d2-p16-d16"
"DELTA2=d2-p16-d16-d12-4u-de+0.6366*p1"


"acqt0=de"
baseopt_echo


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
  (p3 ph1):f2    ; crush eq'm 13C magnetisation
  d13
  p16:gp1
  d16


  (p1 ph1):f1

  DELTA1
  p16:gp2
  d16

  (p3 ph3):f2

;begin purge element
d3
p17:gp3
d16
(center (p2 ph1):f1 (p4 ph1):f2 )
d3
p17:gp3
d16
(center (p2 ph1):f1 (p3 ph6):f2 )
;end purge element

  d0
  (p2 ph5):f1
  d0

  (p3 ph4):f2
  d12 pl12:f2
  p16:gp2
  d16
  4u BLKGRAD
  DELTA2
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2
       F1PH(ip4 & ip29, id0)
       ;F1PH(calph(ph4, +90) & calph(ph29, +90), caldel(d0, +in0))
  4u BLKGRAD
exit 
  

ph1= 0 
ph2= 0 
ph3= 0 2
ph4= 0
ph5= 0 0 0 0 1 1 1 1
ph6= 1 1 3 3
ph29=0
ph31=0 2 0 2 2 0 2 0 


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

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.32

                                          ;preprocessor-flags-start
;SINGLEDWELL: for initial sampling delay of one dwell-time with 
;	    option -DSINGLEDWELL (eda: ZGOPTNS)
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
;	    option -DOFFRES_PRESAT (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: hmqcphpr,v 1.4 2007/04/11 13:34:30 ber Exp $
