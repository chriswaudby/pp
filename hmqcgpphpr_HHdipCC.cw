; For extraction of methyl order parameters by measurement of
; intra-methyl proton-proton dipolar cross-correlated relaxation rates.
; Lit. Tugarinov, Sprangers & Kay, JACS 129, 1743-1750 (2007)
;
; Use DOUBLE_QUANTUM flag to record DQ build-up experiment
; Otherwise, sequence measures decay of 1H SQ magnetisation  

; Delays adjusted for zero first-order phase correction
;
; use baseopt
;_____________________________________________________________
;
;hmqcphpr
;avance-version (07/04/04)
;HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive
;with decoupling during acquisition
;
;solvent suppresion by selective pulse and crusher gradient
;Chris Waudby, April 2013
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

"d11=30m"
"d12=20u"
"d13=4u"

"p2=p1*2"
"p4=p3*2"
"d2=1.8m"		; tau_a
"d3=1m"			; tau_b


"in0=inf2/2"

"d0=in0/2-0.63662*p3-p1"

"DELTA1=2*d2-p16-d16"
"DELTA2=2*d2-p16-d16-d12-de+0.6366*p1+p3"
"DELTA3=d3-p16-d16-larger(p1,p3)"
"acqt0=de"
aqseq 312


1 ze 
2 d11 do:f2

  d12
   "TAU1=vd/4-p3-0.5*p1"
  d12
   "TAU2=vd/4-p19-d16-p3-0.5*p1"

  4u BLKGRAD

3 d12 ;pl9:f1
  d1 ;cw:f1 ph29
  d13 ;do:f1
  d12 pl0:f1 pl2:f2
  4u UNBLKGRAD

  (p11:sp2 ph5):f1    ; SOLVENT SUPPRESSION FLIP-DOWN
  2u 

  (p3 ph6):f2  ; crush equilibrium carbon magnetisation
  d13 pl1:f1
  p16:gp1
  d16 

  p1 ph1
  TAU1
  (p4 ph6):f2
  p19:gp2
  d16
  TAU2
  p2 ph1
  p19:gp2
  d16
  TAU2
  (p4 ph6):f2
  TAU1

#   ifdef DOUBLE_QUANTUM
  p1 ph1
  4u
  p1 ph2
#   else
  p1 ph2
#   endif /*DOUBLE_QUANTUM*/
  
  DELTA1
  p16:gp3
  d16
  (p3 ph3):f2
  p16:gp4
  d16
  DELTA3
  (center (p2 ph6) (p4 ph6):f2 )
  p16:gp4
  d16
  DELTA3
  (p3 ph4):f2
  
  d0
  p2 ph6
  d0

  (p3 ph5):f2
  4u
  p2 ph6

  d12 pl12:f2
  p16:gp3
  d16
  4u BLKGRAD
  DELTA2
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 
	F1QF(ivd)
	F2PH(ip5 & ip29, id0)
  4u BLKGRAD
exit 
  

ph3= 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph4= 1 1 1 1 3 3 3 3 
ph5= 0
ph6= 0
ph29=0

#   ifdef DOUBLE_QUANTUM
ph1= 0 1 2 3
ph2= 0 0 0 0 2 2 2 2 
ph31=0 2 0 2 2 0 2 0 2 0 2 0 0 2 0 2
#   else
ph1= 0 2
ph2= 1 1 3 3
ph31=0 2 0 2 0 2 0 2 2 0 2 0 2 0 2 0
#   endif /*DOUBLE_QUANTUM*/


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;sp2: f1 channel - shaped pulse  90 degree (on H2O)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse 			     [500 usec]
;p19: short homospoil/gradient pulse 		     [50 usec]
;d0 : incremented delay (2D)                         
;d1 : relaxation delay; 1-5 * T1
;d2 : < 1/(4J)CH				     [1.8 msec]
;d3 : delay for purge element - 1/(8J)CH	     [1 msec]
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;vd : variable delay, taken from vd-list
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 16 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence

;use gradient ratio:    gp 1 : gp 2 : gp3 : gp4
;                         18 :  -36 : 18  : 7


;for z-only gradients:
;gpz1: 18%
;gpz2: -36%
;gpz3: 18%
;gpz4: 7%

;use gradient files:
;gpnam1: SINE.10
;gpnam2: SINE.10
;gpnam3: SINE.10
;gpnam4: SINE.10


                                          ;preprocessor-flags-start
;HALFDWELL: for initial sampling delay of half a dwell-time with 
;	    option -DHALFDWELL (eda: ZGOPTNS)
;DOUBLE_QUANTUM: for measurement of build-up of DQ coherence with
;	    option -DDOUBLE_QUANTUM (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: hmqcphpr,v 1.4 2007/04/11 13:34:30 ber Exp $
