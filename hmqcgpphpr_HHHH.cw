; For extraction of methyl order parameters by measurement of
; intra-methyl proton-proton dipolar cross-correlated relaxation rates
; via triple quantum coherence
;
; Lit. Sun, Kay & Tugarinov (2011) J Phys Chem B 115 14878-14884
;
; Use TQ flag to record TQ build-up experiment (ns % 24)
; Use DQ flag to record DQ build-up experiment (ns % 16)
; Otherwise, sequence measures decay of 1H SQ magnetisation (ns % 16)

; can be run as pseudo-2D with -DONE_D (set td2=1)

; Delays adjusted for zero first-order phase correction
; With options for 90,-180 or 180,-360 phase corr.
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

#   ifdef SINGLEDWELL
    "d0=in0-0.63662*p3-p1"
#   else
    "d0=in0/2-0.63662*p3-p1"
#   endif /*SINGLEDWELL*/

"DELTA1=2*d2-p16-d16"
"DELTA2=2*d2-p16-d16-d12-de+0.6366*p1-p3-0.1u"
"DELTA3=d3-p16-d16-larger(p1,p3)"

define delay vdmin
"vdmin=4*(p19+d16+p3+0.5*p1)"

"acqt0=de"

aqseq 312


1 ze 
  vdmin
2 d11 do:f2

  d12
   "TAU1=vd/4-p3-0.5*p1"
  d12
   "TAU2=vd/4-p19-d16-p3-0.5*p1"
  4u UNBLKGRAD


  ; purge before d1
  d12 pl8:f1
  40m cw:f1 ph12  ; 10 kHz purge
  4u do:f1
  p16:gp5
  d16 pl1:f1 
  (p1 ph11)
  p16:gp6
  d16 BLKGRAD

  d12 pl9:f1 fq=cnst21(bf hz):f1
  d1 cw:f1 ph11
  d13 do:f1
  d12 pl2:f2 fq=0:f1
  4u UNBLKGRAD

  ;(p11:sp2 ph11):f1    ; SOLVENT SUPPRESSION FLIP-DOWN
  ;2u 
  (p3 ph11):f2  ; crush equilibrium carbon magnetisation
  d13 pl1:f1
  p16:gp1
  d16 

  p1 ph1
  TAU1
  (p4 ph11):f2
  p19:gp2
  d16
  TAU2
  p2 ph1
  p19:gp2
  d16
  TAU2
  (p4 ph11):f2
  TAU1

#if defined (TQ) || defined (DQ)
  p1 ph2
  0.1u
  p1 ph3
#else
  p1 ph3
#endif /*TQ or DQ */
  
  DELTA1
  p16:gp3
  d16
  (p3 ph4):f2
  p16:gp4
  d16
  DELTA3
  (center (p2 ph11) (p4 ph11):f2 )
  p16:gp4
  d16
  DELTA3

#ifdef ONE_D
  (center (p2 ph11):f1 (p3 ph5 1u p3 ph6):f2 )
#else
  (p3 ph5):f2 
  d0
  p2 ph11
  d0
  (p3 ph6):f2
#endif

  4u
  p2 ph11

  p16:gp3
  d16
  DELTA2
  (p3 ph11):f2
  4u
  (p3 ph7):f2
  d12 pl12:f2 BLKGRAD
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 
	F1QF(ivd)
	F2PH(ip6, id0)
  4u BLKGRAD
exit 
  

ph6= 0
ph11=0
ph12=1

#ifdef TQ
ph1= (12) 0 2 4 6 8 10
ph2= (12) 3 5 7 9 11 1
ph3= 1 1 1 1 1 1 3 3 3 3 3 3
ph4= 0 0 0 0 0 0 0 0 0 0 0 0
     2 2 2 2 2 2 2 2 2 2 2 2
ph5= 1 1 1 1 1 1 3 3 3 3 3 3
ph7= 0 0 0 0 0 0 2 2 2 2 2 2
ph31=0 2 0 2 0 2 0 2 0 2 0 2
     2 0 2 0 2 0 2 0 2 0 2 0
#else
#ifdef DQ
ph1= 0 1 2 3
ph2= 0 1 2 3
ph3= 0 0 0 0 2 2 2 2
ph4= 0 0 0 0 0 0 0 0
     2 2 2 2 2 2 2 2
ph5= 1 1 1 1 3 3 3 3
ph7= 0 0 0 0 0 0 0 0
     2 2 2 2 2 2 2 2
ph31=0 2 0 2 2 0 2 0
     2 0 2 0 0 2 0 2
#else /*SQ*/
ph1= 0 2
ph3= 1 1 3 3
ph4= 0 0 0 0 0 0 0 0
     2 2 2 2 2 2 2 2
ph5= 1 1 1 1 3 3 3 3
ph7= 0 0 0 0 0 0 0 0
     2 2 2 2 2 2 2 2
ph31=0 2 0 2 0 2 0 2
     2 0 2 0 2 0 2 0
#endif
#endif


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl8 : f1 channel - 10 kHz (purge)
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
;cnst21: frequency for off-resonance decoupling (bf Hz)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 24 * n [TQ] / 16 * n [SQ]
;DS: 24 * n [TQ] / 16 * n [SQ]
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
;gpz5: 30%
;gpz6: 37%

;use gradient files:
;gpnam1: SINE.10
;gpnam2: SINE.10
;gpnam3: SINE.10
;gpnam4: SINE.10
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100


                                          ;preprocessor-flags-start
;SINGLEDWELL: for initial sampling delay of half a dwell-time with 
;	    option -DSINGLEDWELL (eda: ZGOPTNS)
;DQ/TQ: for measurement of build-up of DQ / TQ coherence with
;	    option -DDQ / -DTQ (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: hmqcphpr,v 1.4 2007/04/11 13:34:30 ber Exp $
