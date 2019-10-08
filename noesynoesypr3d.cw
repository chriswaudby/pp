;noesynoesypr3d.cw
; with option for 13C/15N decoupling
;avance-version (12/01/11)
;3D sequence with
;   homonuclear correlation via dipolar coupling 
;   dipolar coupling may be due to noe or chemical exchange.
;phase sensitive (t1 and t2)
;with presaturation during relaxation delay and mixing time
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>


"d11=30m"
"d12=20u"
"d13=4u"


"in0=inf1"
"in10=inf2"

"d0=in0/2-p1*4/3.1416-1u"
"d10=in10/2-p1*4/3.1416-1u"

"acqt0=-p1*2/3.1416"

aqseq 312


1 d11 ze
#ifdef LABEL_CN
  d11 pl12:f2 pl26:f3
2 d11 do:f2 do:f3
#else
  d11
2 d11
#endif
  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  p1 ph1
#ifdef LABEL_CN
  d0 cpd2:f2 cpd3:f3
  1u do:f2 do:f3
#else
  d0
  1u
#endif
  p1 ph2
  d12 pl9:f1
  d8 cw:f1
  d13 do:f1
  d12 pl1:f1
  p1 ph3
#ifdef LABEL_CN
  d10 cpd2:f2 cpd3:f3
  1u do:f2 do:f3
#else
  d10
  1u
#endif
  p1 ph4
  d12 pl9:f1
  d20 cw:f1
  d13 do:f1
  d12 pl1:f1
  p1 ph5
#ifdef LABEL_CN
  go=2 ph31 cpd2:f2 cpd3:f3
  d11 do:f2 do:f3 mc #0 to 2 
#else
  go=2 ph31
  f11 mc #0 to 2
#endif
     F1PH(calph(ph1, +90) & calph(ph29, +90), caldel(d0, +in0)) 
     F2PH(calph(ph3, +90), caldel(d10, +in10))
#ifdef LABEL_CN
  d11 do:f2 do:f3
#endif
exit 
  

ph1=0 2
ph2=0
ph3=0
ph4=0 0 2 2
ph5=0
ph29=0
ph31=0 2 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;p1 : f1 channel -  90 degree high power pulse
;d0 : incremented delay (F1 in 3D)
;d1 : relaxation delay; 1-5 * T1
;d8: first mixing time
;d10: incremented delay (F2 in 3D)
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d20: second mixing time
;inf1: 1/SW(H) = 2 * DW(H)
;inf2: 1/SW(H) = 2 * DW(H)
;in0: 1/(1 * SW(H)) = 2 * DW(H)
;nd0: 1
;in10: 1/(1 * SW(H)) = 2 * DW(H)
;nd10: 1
;ns: 4 * n
;ds: 32
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI (or TPPI) in F1
;FnMODE: States-TPPI (or TPPI) in F2

;for older datasets use AQORDER : 3 - 1 - 2


;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1

;PHC0(F2): 90
;PHC1(F2): -180
;FCOR(F2): 1



;$Id: noesynoesypr3d,v 1.5 2012/01/31 17:49:28 ber Exp $

