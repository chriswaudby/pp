;noesyhsqcgpsm3d.2
;avance-version (12/01/11)
;NOESY-HSQC
;3D sequence with
;   homonuclear correlation via dipolar coupling 
;   dipolar coupling may be due to noe or chemical exchange.
;   H-1/X correlation via double inept transfer
;   simultaneous evolution of C-13 and N-15 chemical shift in t2
;phase sensitive (t1)
;phase sensitive (t2)
;with decoupling during acquisition
;using shaped pulses for inversion on f2 - channel
;(use parameterset NOESYHSQCGPSM3D.2)
;
;S.M. Pascal, D.R. Muhandiram, T. Yamazaki, J.D. Forman-Kay & L.E. Kay, 
;   J. Magn. Reson. B103, 197 - 201 (1994)
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
"p4=p3*2"
"p22=p21*2"
"d4=1s/(cnst2*4)"
"d11=30m"
"d12=20u"
"d26=1s/(cnst4*4)"
"d29=3u"

"p16=1m"
"p19=3m"
"p29=4m"


"d0=3u"
"d10=3u"

"in0=inf1/2"
"in10=inf2/2"

"in29=in30/2-in10"


"DELTA=d1-d11*4-6m"
"DELTA1=d4-p16-d16-p8/2"
"DELTA2=d26-d4"
"DELTA3=d4-p16-d16-p8/2-p3*2-p21*2-d12-10u"
"DELTA4=d0*2+larger(p8,p22)"
"DELTA5=p2*2+p3*2+p4+d10*4+d29*2"
"DELTA6=d10*2+p2"

"TAU=d8-p19*2-d16*2-p3"
"TAU2=(p8-p2)/2"
"TAU3=(p8-p22)/2"


"spoff13=0"


aqseq 321


1 ze
  d11 pl12:f2 pl16:f3
2 d11 do:f2 do:f3
3 d12 pl10:f1
  p17 ph1
  p17*2 ph2

  d1 pl1:f1 pl0:f2 pl3:f3
  50u UNBLKGRAD
  p16:gp7
  d16

  (p1 ph3)
  DELTA4
  (p2 ph8)
  d0
  (center (p8:sp13 ph1):f2 (p22 ph1):f3 )
  d0
  (p1 ph1)

  TAU
  p19:gp1
  d16 pl2:f2
  (p3 ph1):f2
  4u
  p16:gp2
  d16

  (p1 ph1)
  4u
  p16:gp3
  d16 pl0:f2
  DELTA1
  (lalign
   (DELTA2 TAU2 p2 ph1 TAU2)
   (p8:sp13 ph1):f2
   (DELTA2 TAU3 p22 ph1 TAU3):f3
  )
  DELTA2
  DELTA1
  4u
  p16:gp3
  d16
  (p1 ph2)

  4u
  p29:gp4
  d16 pl2:f2

  (p21 ph5):f3
  d29
  (p3 ph4):f2 
  d10
  (p2 ph7)
  d10
  (p4 ph9):f2 
  DELTA6
  (p3 ph1):f2
  d29
  (p22 ph9):f3
  DELTA5
  (p21 ph1):f3 

  4u
  p16:gp5
  d16 

  (p1 ph6)
  4u
  p16:gp6
  d16 pl0:f2
  DELTA1
  (lalign
   (DELTA2 TAU2 p2 ph1 TAU2):f1
   (p8:sp13 ph1):f2
   (DELTA2 TAU3 p22 ph1 TAU3):f3
  )
  DELTA2
  DELTA3
  4u
  p16:gp6
  d16 pl2:f2
  (p3 ph1 3u p3 ph7):f2
  (p21 ph1 3u p21 ph7):f3
  d12 pl12:f2 pl16:f3
  4u BLKGRAD
  (p1 ph6)

  go=2 ph31 cpd2:f2 cpd3:f3
  d11 do:f2 do:f3 mc #0 to 2 
     F1PH(calph(ph3, +90) & calph(ph8, +90), caldel(d0, +in0)) 
     F2PH(calph(ph4, +90) & calph(ph5, +90), caldel(d10, +in10) & caldel(d29, +in29))
exit
   

ph1=0 
ph2=1
;ph3=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph3=0 0 0 0 2 2 2 2
ph4=0 2
ph5=2 0
;ph6=0 0 1 1 2 2 3 3
ph6=0 0 1 1 0 0 1 1 2 2 3 3 2 2 3 3
;ph7=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph7=0 0 0 0 2 2 2 2
;ph8=1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3
ph8=1 1 1 1 3 3 3 3
ph9=0 0 1 1
ph31=0 2 1 3 2 0 3 1 2 0 3 1 0 2 1 3


;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl10: f1 channel - power level for TOCSY-spinlock (trim pulse)
;pl12: f2 channel - power level for CPD/BB decoupling
;pl16: f3 channel - power level for CPD/BB decoupling
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse                       [1 msec]
;p17: f1 channel - trim pulse                        [4 msec]
;p19: gradient pulse 2                               [3 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p29: gradient pulse 3                               [4 msec]
;d0 : incremented delay (F1 in 3D)                   [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J(CH))
;d8 : mixing time
;d10: incremented delay (F2 in 3D)                   [3 usec]
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d26: 1/(4J(NH))
;d29: decremented delay (F2 in 3D)                   [3 usec]
;cnst2: = J(CH)
;cnst4: = J(NH)
;inf1: 1/SW(H) = 2 * DW(H)
;inf2: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(H)) =  DW(H)
;nd0: 2
;in10: 1/(2 * SW(C)) = DW(C)
;nd10: 2
;in29: in30/2 - in10
;in30: 1/SW(N) = 2 * DW(N)
;ns: 8 * n
;ds: >= 16
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI (or TPPI) in F1
;FnMODE: States-TPPI (or TPPI) in F2
;cpd2: decoupling according to sequence defined by cpdprg2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    gp 1 : gp 2 : gp 3 : gp 4 : gp 5 : gp 6 : gp 7
;                         30 :   40 :   16 :   60 :  -36 :   16 :   50

;for z-only gradients:
;gpz1: 30%
;gpz2: 40%
;gpz3: 16%
;gpz4: 60%
;gpz5: -36%
;gpz6: 16%
;gpz7: 50%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam7: SMSQ10.100



;$Id: noesyhsqcgpsm3d.2,v 1.9 2012/01/31 17:49:28 ber Exp $
