;Recoded 1H indirect shift evolution to allow for wide SW without 180deg pulse
;15N coherence-order-selection gradient implemented as bipolar pair
;Delays adjusted for zero first-order phase correction
;
;dipsihsqcf3gpsi3d
;avance-version (07/04/04)
;TOCSY-HSQC
;3D sequence with
;   homonuclear Hartman-Hahn transfer using DIPSI2 sequence
;      for mixing
;   H-1/X correlation via double inept transfer
;      using sensitivity improvement
;phase sensitive (t1)
;phase sensitive using Echo/Antiecho-TPPI gradient selection (t2)
;using trim pulses in inept transfer
;using f3 - channel
;(use parameterset DIPSIHSQCF3GPSI3D)
;
;A.G. Palmer III, J. Cavanagh, P.E. Wright & M. Rance, J. Magn.
;   Reson. 93, 151-170 (1991)
;L.E. Kay, P. Keifer & T. Saarinen, J. Am. Chem. Soc. 114,
;   10663-5 (1992)
;J. Schleucher, M. Schwendinger, M. Sattler, P. Schmidt, O. Schedletzky,
;   S.J. Glaser, O.W. Sorensen & C. Griesinger, J. Biomol. NMR 4,
;   301-306 (1994)
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"
"d13=4u"
"d26=1s/(cnst4*4)"



"in0=inf1"
"in10=inf2/2"

# ifdef RESUME
    "d0=in0/2-1.27324*p1+l31*in0"
# else
    "d0=in0/2-1.27324*p1"
# endif /*RESUME*/

"d10=3u"

"DELTA1=d13+p16+d16+4u+de-0.63662*p1"

# ifdef LABEL_CN
   "DELTA=d10*2+larger(p2,p14)-4u"
# else
   "DELTA=d10*2+p2-4u"
# endif /*LABEL_CN*/


"FACTOR1=(d9/(p6*115.112))/2+0.5"
"l1=FACTOR1*2"


aqseq 321


1 ze
  d11 pl16:f3
2 d11 do:f3
3 d12 pl9:f1
  d1 cw:f1 ph29
  4u do:f1
  d12 pl1:f1 pl3:f3

# ifdef LABEL_CN
  
  d12 pl0:f2

  if "d0 < p14" goto 4
  10u
  (center (p1 ph8 d0 p1 ph10):f1 (p14:sp3 ph1):f2 (p22 ph1):f3 )
  goto 6

4 10u
  if "d0 < p22" goto 5
  (center (p1 ph8 d0 p1 ph10):f1 (p22 ph1):f3 )
  goto 6

5 (p1 ph8 d0 p1 ph10):f1

# else
  
  d12
  10u
  if "d0 < p22"
  {
  (p1 ph8 d0 p1 ph10):f1
  }
  else
  {
  (center (p1 ph8 d0 p1 ph10):f1 (p22 ph1):f3)
  }

# endif /*LABEL_CN*/

6 10u
  d20 pl10:f1

						;begin DIPSI2
7 p6*3.556 ph23
  p6*4.556 ph25
  p6*3.222 ph23
  p6*3.167 ph25
  p6*0.333 ph23
  p6*2.722 ph25
  p6*4.167 ph23
  p6*2.944 ph25
  p6*4.111 ph23
  
  p6*3.556 ph25
  p6*4.556 ph23
  p6*3.222 ph25
  p6*3.167 ph23
  p6*0.333 ph25
  p6*2.722 ph23
  p6*4.167 ph25
  p6*2.944 ph23
  p6*4.111 ph25

  p6*3.556 ph25
  p6*4.556 ph23
  p6*3.222 ph25
  p6*3.167 ph23
  p6*0.333 ph25
  p6*2.722 ph23
  p6*4.167 ph25
  p6*2.944 ph23
  p6*4.111 ph25

  p6*3.556 ph23
  p6*4.556 ph25
  p6*3.222 ph23
  p6*3.167 ph25
  p6*0.333 ph23
  p6*2.722 ph25
  p6*4.167 ph23
  p6*2.944 ph25
  p6*4.111 ph23
  lo to 7 times l1
						;end DIPSI2

  d21 pl1:f1
  (p1 ph11)

  d26
  (center (p2 ph1) (p22 ph6):f3 )
  d26 UNBLKGRAD
  p28 ph1
  d13
  (p1 ph2) 
  3u
  p16:gp1
  d16
  (p21 ph3):f3
  d10 

# ifdef LABEL_CN
   (center (p2 ph7):f1 (p14:sp3 ph1):f2 )
# else
   (p2 ph7):f1
# endif /*LABEL_CN*/

  d10
  p16:gp2*EA
  d16
  (p22 ph4):f3
  4u
  p16:gp2*EA*-1
  d16
  DELTA
  
  (center (p1 ph1) (p21 ph4):f3 )
  d24
  (center (p2 ph1) (p22 ph1):f3 )
  d24
  (center (p1 ph2) (p21 ph5):f3 )
  d26
  (center (p2 ph1) (p22 ph1):f3 )
  d26
  (p1 ph1)
  DELTA1
  (p2 ph1)
  d13
  p16:gp3
  d16 pl16:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
     F1PH(rd10 & rp3 & rp6 & rp31 & ip8 & ip9 & ip29, id0) 
     F2EA(igrad EA & ip5*2, id10 & ip3*2 & ip6*2 & ip31*2)
exit
   

ph1=0 
ph2=1
ph3=0 2
ph4=0 0 2 2
ph5=1 1 3 3
ph6=0
ph7=0 0 2 2
ph8=0 0 0 0 2 2 2 2
ph10=2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
ph11=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph23=1
ph25=3
ph29=0
ph31=0 2 2 0 2 0 0 2 2 0 0 2 0 2 2 0
     2 0 0 2 0 2 2 0 0 2 2 0 2 0 0 2
  

;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl10: f1 channel - power level for TOCSY-spinlock
;pl16: f3 channel - power level for CPD/BB decoupling
;sp3 : f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p6 : f1 channel -  90 degree low power pulse
;p14: f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse                         [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p28: f1 channel - trim pulse                          [1 msec]
;d0 : incremented delay (F1 in 3D)                     [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d9 : TOCSY mixing time
;d10: incremented delay (F2 in 3D)                     [3 usec]
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d13: short delay                                      [4 usec]
;d16: delay for homospoil/gradient recovery
;d20: first z-filter delay                           [10 usec]
;d21: second z-filter delay                          [10 usec]
;d24: 1/(4J)YH for YH
;     1/(8J)YH for all multiplicities
;d26: 1/(4J(YH))
;cnst4: = J(YH)
;l1: loop for DIPSI cycle: ((p6*115.112) * l1) = mixing time
;l31: completed complex pairs in F1 for resuming after interruption
;inf1: 1/SW(H) = 2 * DW(H)
;inf2: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(H)) = DW(H)
;nd0: 2
;in10: 1/(2 * SW(X)) = DW(X)
;nd10: 2
;NS: 8 * n
;DS: >= 16
;td1: number of experiments
;td2: number of experiments in F2
;FnMODE: States-TPPI (or TPPI) in F1
;FnMODE: echo-antiecho in F2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    gp 1 : gp 2 : gp 3
;                         50 :   80 : 16.21

;for z-only gradients:
;gpz1: 50%
;gpz2: 80%
;gpz3: 16.21%

;use gradient files:   
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100


;set pl9 to 120dB when presaturation is not required
;   use 70 - 80dB to reduce radiation damping


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: dipsihsqcf3gpsi3d,v 1.4 2007/04/11 13:34:29 ber Exp $
