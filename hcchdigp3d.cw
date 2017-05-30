;hcchdigp3d
;avance-version (14/02/11)
;HCCH-TOCSY
;3D sequence with
;   inverse correlation using multiple inept transfer and
;      C-C DIPSI3 spinlock
;
;      F1(H,t1) -> F2(C,t2) -> F2(C') -> F1(H',t3)
;
;off resonance C=O pulse using shaped pulse
;phase sensitive (t1)
;phase sensitive (t2)
;spinlock during z-filter
;(use parameterset HCCHDIGP3D)
;
;(L.E. Kay, G.Y. Xu, A.U. Singer, D.R. Muhandiram & J. D. Forman-Kay
;   J. Magn. Reson. B 101, 333 - 337 (1993))
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
"d11=30m"
"d12=20u"

"d4=1.6m"			;tau a
"d21=1.1m"			;tau c
"d23=475u"			;tau b


"p16=500u"
"p19=2m"
"p29=300u"
"p30=5m"
"p31=4.4m"

"d16=150u"


"d0=3u"
"d10=3u"

"in0=inf1/2"
"in10=inf2/2"


"DELTA1=d4-p16-d16-4u"
"DELTA2=d4-p16-d16-4u+d0*2+p4"
"DELTA3=d23-p29-d16"
"DELTA4=p22+p2+d10*2+4u"
"DELTA5=d21-p16-d16-4u"
"DELTA6=d4-p16-d16-p3*2-15u"


"d31=(p9*54.33*4)*l1"


"spoff5=bf2*(cnst21/1000000)-o2"


"acqt0=-p1*2/PI"


aqseq 312


1 d11 ze
  d31 pl12:f2  pl3:f3
2 d11 do:f2
3 d1
  50u UNBLKGRAD
  d12 pl1:f1

  (p1 ph3)
  4u
  p16:gp1
  d16
  DELTA1 pl2:f2
  d0
  (p4 ph1):f2
  d0
  (p2 ph1)
  4u
  p16:gp1
  d16
  DELTA2 pl3:f3
  (p1 ph2)

  p19:gp3
  d16

  (p3 ph4):f2
  d10
  (p22 ph1):f3
  4u
  p29:gp2
  d16
  DELTA3 pl0:f2
  (p14:sp5 ph1):f2
  4u
  p29:gp2
  d16
  DELTA3 pl2:f2
  p2 ph1
  d10
  (p4 ph1):f2
  DELTA4
  p29:gp2
  d16
  DELTA3
  (p14:sp5 ph1):f2
  4u
  p29:gp2
  d16
  DELTA3 pl2:f2
  (p3 ph2):f2
  4u
  d12 pl15:f2
						;begin DIPSI3
9 (p9*2.722 ph7):f2
  (p9*4.389 ph9):f2
  (p9*2.778 ph7):f2
  (p9*3.056 ph9):f2
  (p9*0.333 ph7):f2
  (p9*2.556 ph9):f2
  (p9*4.000 ph7):f2
  (p9*2.722 ph9):f2
  (p9*4.111 ph7):f2
  (p9*3.778 ph9):f2
  (p9*3.889 ph7):f2
  (p9*2.889 ph9):f2
  (p9*3.000 ph7):f2
  (p9*0.333 ph9):f2
  (p9*2.500 ph7):f2
  (p9*4.050 ph9):f2
  (p9*2.830 ph7):f2
  (p9*4.389 ph9):f2
  (p9*2.722 ph9):f2
  (p9*4.389 ph7):f2
  (p9*2.778 ph9):f2
  (p9*3.056 ph7):f2
  (p9*0.333 ph9):f2
  (p9*2.556 ph7):f2
  (p9*4.000 ph9):f2
  (p9*2.722 ph7):f2
  (p9*4.111 ph9):f2
  (p9*3.778 ph7):f2
  (p9*3.889 ph9):f2
  (p9*2.889 ph7):f2
  (p9*3.000 ph9):f2
  (p9*0.333 ph7):f2
  (p9*2.500 ph9):f2
  (p9*4.050 ph7):f2
  (p9*2.830 ph9):f2
  (p9*4.389 ph7):f2
  (p9*2.722 ph9):f2
  (p9*4.389 ph7):f2
  (p9*2.778 ph9):f2
  (p9*3.056 ph7):f2
  (p9*0.333 ph9):f2
  (p9*2.556 ph7):f2
  (p9*4.000 ph9):f2
  (p9*2.722 ph7):f2
  (p9*4.111 ph9):f2
  (p9*3.778 ph7):f2
  (p9*3.889 ph9):f2
  (p9*2.889 ph7):f2
  (p9*3.000 ph9):f2
  (p9*0.333 ph7):f2
  (p9*2.500 ph9):f2
  (p9*4.050 ph7):f2
  (p9*2.830 ph9):f2
  (p9*4.389 ph7):f2
  (p9*2.722 ph7):f2
  (p9*4.389 ph9):f2
  (p9*2.778 ph7):f2
  (p9*3.056 ph9):f2
  (p9*0.333 ph7):f2
  (p9*2.556 ph9):f2
  (p9*4.000 ph7):f2
  (p9*2.722 ph9):f2
  (p9*4.111 ph7):f2
  (p9*3.778 ph9):f2
  (p9*3.889 ph7):f2
  (p9*2.889 ph9):f2
  (p9*3.000 ph7):f2
  (p9*0.333 ph9):f2
  (p9*2.500 ph7):f2
  (p9*4.050 ph9):f2
  (p9*2.830 ph7):f2
  (p9*4.389 ph9):f2
  lo to 9 times l1
						;end DIPSI3
  d12 pl10:f1
  (p17 ph1)
  (p17*2 ph2)
  4u
  p30:gp4
  d16 pl1:f1
  (p1 ph1)
  4u
  p31:gp4
  d16 pl2:f2

  (p3 ph2):f2
  4u
  p16:gp1
  d16
  DELTA5
  (center (p2 ph1) (p4 ph1):f2 )
  4u
  p16:gp1
  d16
  DELTA5
  (p3 ph1):f2

  (p1 ph1)
  4u
  p16:gp1
  d16
  DELTA1
  (center (p2 ph1) (p4 ph1):f2 )
  4u
  p16:gp1
  d16
  DELTA6 
  4u BLKGRAD
  (p3 ph1 3u p3 ph5):f2
  4u pl12:f2 
  (p1 ph1)
  go=2 ph31 cpd2:f2
  d11 do:f2 mc #0 to 2
     F1PH(calph(ph3, +90), caldel(d0, +in0) & calph(ph3, -90))
     F2PH(calph(ph4, +90), caldel(d10, +in10))
exit


ph1=0
ph2=1
ph3=0 2
ph4=0 0 2 2
ph5=0 0 0 0 2 2 2 2
ph7=1
ph9=3
ph31=0 2 2 0


;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl10: f1 channel - power level for TOCSY-spinlock (trim pulse)
;pl12: f2 channel - power level for CPD/BB decoupling
;pl15: f2 channel - power level for TOCSY-spinlock
;sp5: f2 channel - shaped pulse 180 degree   (C=O off resonance)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p9 : f2 channel -  90 degree low power pulse
;p14: f2 channel - 180 degree shaped pulse
;p16: homospoil/gradient pulse                         [500 usec]
;p17: f1 channel - trim pulse                          [2.5 msec]
;p19: gradient pulse 2                                 [2 msec]
;p22: f3 channel - 180 degree high power pulse
;p29: gradient pulse 3                                 [300 usec]
;p30: gradient pulse 4                                 [5 msec]
;p31: gradient pulse 5                                 [4.4 msec]
;d0 : incremented delay (F1 in 3D)                     [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J(CH)) - tau a                               [1.6 msec]
;d10: incremented delay (F2 in 3D)                     [3 usec]
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d16: delay for homospoil/gradient recovery
;d21: 1/(6J'(CH)) - tau c                              [1.1 msec]
;d23: tau b                                            [475 usec]
;d31: length of DIPSI-3 cycle as executed = (p9*54.33*4)*l1
;cnst21: CO chemical shift (offset, in ppm)
;cnst23: Caliphatic chemical shift (offset, in ppm)
;o2p: Caliphatic chemical shift (cnst23)
;l1: loop for DIPSI cycle:
;       mixing time = ((p9*54.33*4) * l1)              [12 msec]
;inf1: 1/SW(Hali) = 2 * DW(Hali)
;inf2: 1/SW(C) = 2 * DW(C)
;in0: 1/(2 * SW(Hali)) = DW(Hali)
;nd0: 2
;in10: 1/(2 * SW(C)) = DW(C)
;nd10: 2
;ns: 8 * n
;ds: 32
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI in F1 (not suitable for TPPI)
;FnMODE: States-TPPI (or TPPI) in F2
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    gp 1 : gp 2 : gp 3 : gp 4
;                         16 :   16 :   30 :   60

;for z-only gradients:
;gpz1: 16%
;gpz2: 16%
;gpz3: 30%
;gpz4: 60%

;use gradient files:
;gpnam1: SMSQ10.50
;gpnam2: SMSQ10.50
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100



;Processing

;SR(F1): 1/4 SWH(F1)



;$Id: hcchdigp3d,v 1.19.2.1 2014/02/11 11:32:03 ber Exp $
