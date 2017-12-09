;Clean-up gradient pairs added in SE block, bipolar gradients during t1
;Delays adjusted for zero first-order phase correction
;
;hsqct2etf3gpsi3d
;avance-version (07/04/04)
;3D H-1/X correlation via double inept transfer
;   using sensitivity improvement
;for measuring N-15 T2 relaxation times
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;using f3 - channel
;using flip-back pulse
;as pseudo3D
;(use parameterset HSQCT2ETF3GPSI3D)
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>



"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"
"d24=1s/(cnst4*cnst11)"
"d25=1s/(cnst4*cnst12)"
"d26=1s/(cnst4*4)"

"d21=450u"

"d31=(p30*16+d21*32)"

		
"d10=10u"

"in10=inf2/4"

"DELTA2=p16+d16+8u+d12+de-0.6366*p1"
"DELTA3=d21-p2/2"
"DELTA4=d25-p16-d16"
"DELTA5=d24-p16-d16-4u"
"DELTA6=d26-p16-d16-4u"

#   ifdef LABEL_CN
"DELTA1=d25-p16-d16-larger(p2,p8)-d10*2-p21*4/3.1415"
"spoff13=bf2*((cnst21+cnst22)/2000000)-o2"
#   else
"DELTA1=d25-p16-d16-p2-d10*2-p21*4/3.1415"
#   endif /*LABEL_CN*/


"spoff1=0"


aqseq 312


1 ze
  d11 pl16:f3 st0
2 6m do:f3 
3 3m
4 d1

  (p1 ph1)
  d26 pl3:f3
  (center (p2 ph1) (p22 ph1):f3 )
  d26 UNBLKGRAD
  (p1 ph2)

  4u pl0:f1
  (p11:sp11 ph1:r):f1		; flipback(+x), +y -> +z
  4u
  p16:gp1
  d16 pl1:f1

  (p21 ph3):f3
  d25 
  (center (p2 ph1) (p22 ph6):f3 )	; water -> -z
  d25 pl23:f3

6 d21
  (p30 ph6 d21*2 p30 ph6 d21*2 p30 ph6 d21*2 p30 ph6):f3
  DELTA3
  (p2 ph1)				; water -> +z
  DELTA3
  (p30 ph6 d21*2 p30 ph6 d21*2 p30 ph6 d21*2 p30 ph6):f3
  d21*2
  (p30 ph6 d21*2 p30 ph6 d21*2 p30 ph6 d21*2 p30 ph6):f3
  DELTA3
  (p2 ph8)				; water -> -z
  DELTA3
  (p30 ph6 d21*2 p30 ph6 d21*2 p30 ph6 d21*2 p30 ph6):f3
  d21
  lo to 6 times c

  p16:gp2*-1*EA
  d16
  DELTA4 pl3:f3
  (p22 ph6):f3

  2u
  d10 gron0
  d10 gron0*-1
  8u groff 

#   ifdef LABEL_CN
  (center (p2 ph7) (p8:sp13 ph1):f2 )	; water -> +z
#   else
  (p2 ph7)
#   endif /*LABEL_CN*/

  2u
  d10 gron0
  d10 gron0*-1
  8u groff 

  p16:gp2*EA
  d16
  DELTA1

  (center (p1 ph1) (p21 ph4):f3 )	; water -> -y
  4u
  p16:gp4
  d16
  DELTA5
  (center (p2 ph1) (p22 ph1):f3 )	; water -> +y
  DELTA5
  4u
  p16:gp4
  d16
  (center (p1 ph2) (p21 ph5):f3 )	; water -- +y
  4u
  p16:gp5
  d16
  DELTA6
  (center (p2 ph1) (p22 ph1):f3 )	; water -> -y
  DELTA6
  4u
  p16:gp5
  d16
  (p1 ph1)				; water -> -z
  DELTA2
  (p2 ph1)				; water -> +z
  4u
  p16:gp3
  d16
  d12 pl16:f3
  4u  BLKGRAD
  goscnp ph31 cpd3:f3

  3m do:f3
  3m st ivc
  lo to 3 times nbl

  3m ipp3 ipp4 ipp5 ipp6 ipp7 ipp31
  lo to 4 times ns

  d1 mc #0 to 4
     F1QF()
     F2EA(igrad EA & ip5*2 & rpp3 rpp4 rpp5 rpp6 rpp7 rpp31, id10 & ip3*2 & ip6*2 & ip31*2)
  d31
exit
   

ph0=0 
ph1=0 
ph2=1
ph3=0 2 
ph4=0 0 2 2
ph5=3 3 1 1
ph6=0 0 0 0 2 2 2 2
ph7=0 0 2 2
ph8=2
ph31=0 2 2 0
  

;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl16: f3 channel - power level for CPD/BB decoupling
;pl23: f3 channel - power level for TOCSY-spinlock
;sp1 : f1 channel - shaped pulse  90 degree
;sp13: f2 channel - shaped pulse 180 degree  (Ca and C=O, adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p30: f3 channel - 180 degree pulse at pl23
;d1 : relaxation delay; 1-5 * T1
;d10 : incremented delay                             [3 usec]
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21: echo delay                                     [450 usec]
;d24: 1/(4J)YH for YH
;     1/(8J)YH for all multiplicities
;d25: 1/(4J)YH for YH
;     1/(8J)YH for all multiplicities
;d26: 1/(4J(YH))
;d31: length of single cpmg loop
;cnst4: = J(YH)
;cnst11: for multiplicity selection = 4 for NH, 8 for all multiplicities
;cnst12: for multiplicity selection = 4 for NH, 8 for all multiplicities
;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)
;vc : variable loop-coounter for T2 delay, taken from vc-list
;inf2: 1/SW(X) = 2 * DW(X)
;in10: 1/(2 * SW(X)) = DW(X)
;nd10: 2
;NS: 2 * n
;DS: >= 16
;td1: number of delays in vc-list
;td2: number of experiments in F2
;NBL: = td1
;FnMODE: QF in F1
;FnMODE: echo-antiecho in F2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz0: 1-2%
;gpz1: 30%
;gpz2: 80%
;gpz3: 16.2%
;gpz4: 7%
;gpz5: -5%

;use gradient files:   
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam5: SINE.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: hsqct2etf3gpsi3d,v 1.5 2007/04/11 13:34:30 ber Exp $
