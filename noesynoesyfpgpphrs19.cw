;noesynoesyfpgpphrs19.cw
;avance-version (12/01/11)
;2D homonuclear correlation via dipolar coupling 
;dipolar coupling may be due to noe or chemical exchange.
;phase sensitive
;water suppression using 3-9-19 pulse sequence with gradients
;using flip-back pulse
;with radiation damping suppression using gradients in t1
;
;M. Piotto, V. Saudek & V. Sklenar, J. Biomol. NMR 2, 
;   661 - 666 (1992)
;V. Sklenar, M. Piotto, R. Leppik & V. Saudek, J. Magn. Reson. A102,
;   241 -245 (1993)
;G. Lippens, C. Dhalluin & J.-M. Wieruszeski, J. Biomol. NMR 5,
;   327-331 (1995)
;V. Sklenar, J. Magn. Reson. A114, 132-135 (1995)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"d12=20u"


"in0=inf1/2"
"in10=inf2/2"

"d0=in0/2-p1*2/3.1416-1u"
"d10=in10/2-p1*2/3.1416-1u"

"DELTA1=d8"
"DELTA2=d8-p16-d16-p11-d12-4u"

"TAU=p1*2/3.1416-(p0-p27)*0.231-de+46u"

aqseq 312
"acqt0=0"
baseopt_echo


1 ze
2 d1
3 d12 pl1:f1


  50u UNBLKGRAD
  ; t1 evolution
  p1 ph1
  d0 gron0
  d0 gron0*-1
  2u groff
  p1 ph2

  ; first mixing time
  DELTA1

  ; t2 evolution
  p1 ph11
  d10 gron0
  d10 gron0*-1
  2u groff
  p1 ph2

  ; second mixing time
  DELTA2

  ; final readout
  p16:gp1
  d16 pl0:f1
  (p11:sp1 ph3:r):f1
  4u
  d12 pl1:f1
  p1 ph4
  50u pl18:f1
  p16:gp2
  d16
  p27*0.231 ph5
  d19*2
  p27*0.692 ph5
  d19*2
  p27*1.462 ph5
  d19*2
  p27*1.462 ph6
  d19*2
  p27*0.692 ph6
  d19*2
  p0*0.231 ph6
  TAU
  p16:gp2
  d16
  4u BLKGRAD
  go=2 ph31
  d1 mc #0 to 2
       F1PH(calph(ph1, +90), caldel(d0, +in0))
       F2PH(calph(ph11, +90), caldel(d10, +in0))
exit


ph1= 0 2 
ph11=0 0 2 2
ph2= 0
ph3= 2 2 2 2 0 0 0 0 
ph4= 0 0 0 0 2 2 2 2 
ph5= 0
ph6= 2
ph31=0 2 2 0 2 0 0 2 

;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl18: f1 channel - power level for 3-9-19-pulse (watergate)
;sp1: f1 channel - shaped pulse  90 degree
;p0 : f1 channel -  90 degree pulse at pl18
;                      use for fine adjustment
;p1 : f1 channel -  90 degree high power pulse
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse
;p27: f1 channel -  90 degree pulse at pl18
;d0 : incremented delay (2D, min >= 6usec)
;d1 : relaxation delay; 1-5 * T1
;d8 : mixing time
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d19: delay for binomial water suppression
;     d19 = (1/(2*d)), d = distance of next null (in Hz)
;inf1: 1/SW = 2 * DW
;in0: 1/(2 * SW) = DW
;nd0: 1
;ns: 4 * n
;ds: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ


;use gradient ratio:    gp 0 : gp 1 : gp 2
;                          2 :   50 :   30

;for z-only gradients:
;gpz0: 2%
;gpz1: 50%
;gpz2: 30%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1



;$Id: noesyfpgpphrs19,v 1.12 2012/01/31 17:49:27 ber Exp $
