;hisqcctetgpsp
; 13C in-phase HSQC
; starting with equilibrium 13Cz magnetisation
;
;avance-version (15/02/27)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;constant time version
;using trim pulses in inept transfer
;using shaped pulses for inversion on f2 - channel
;
;(G.W. Vuister & A. Bax, J. Magn. Reson. 98, 428-435 (1992))
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


"p2=p1*2"
"p22=p21*2"
"d3=1s/(cnst2*6)"
"d4=1s/(cnst2*4)"
"d11=30m"
"d12=20u"


"d0=inf1/2"
"in0=inf1"

"TAU=d3-4u-0.5*larger(p14,p2)"
"TAU1=d3-p16-d16-0.5*larger(p14,p2)"

;"td1=tdmax(td1,d20*2,in20)"


"DELTA1=d4-larger(p2,p8)/2-p16-de-8u"
"DELTA2=d4-larger(p2,p8)/2-4u-p16-d16"
;"DELTA3=d3-larger(p2,p8)/2-4u-p19-d16"
"DELTA3=d3-larger(p2,p14)/2-4u-p19-d16"
"DELTA4=d4-larger(p2,p8)/2-p1*2/PI-4u-p16-d16"


"spoff3=0"
"spoff5=bf2*(cnst21/1000000)-o2"
"spoff13=0"


"acqt0=0"
baseopt_echo


1 ze 
  d11 pl12:f2 pl3:f3
2 d1 do:f2

  4u UNBLKGRAD
  4u pl1:f1 pl2:f2

#ifdef C13START
  ; purge equilibrium 1H magnetisation
  (p1 ph1)
  4u
  p16:gp0
  d16
#else
  ; purge equilibrium 13C magnetisation
  (p3 ph1):f2
  4u
  p16:gp0
  d16

  ; first INEPT
  (p1 ph1)
  4u
  p16:gp4
  d16
  DELTA2 pl0:f2
  (center (p2 ph1) (p8:sp13 ph6):f2 )
  DELTA2 pl2:f2
  p16:gp4
  d16
  (p1 ph2)
  ; zz filter
  4u
  p16:gp3*-0.8
  d16
  ; second INEPT
  (p3 ph1):f2
  4u
  p19:gp4
  d16
  DELTA3 pl0:f2
  ;(center (p2 ph1) (p8:sp13 ph6):f2 )
  (center (p2 ph1) (p14:sp3 ph6):f2 )
  DELTA3 pl2:f2
  p19:gp4
  d16
  4u
  (p3 ph2):f2
#endif /*C13START*/

  ; purge water
  4u
  p16:gp0*1.1
  d16 pl6:f1
  (2mp ph1):f1
  (3mp ph2):f1

#ifdef ZFILTER
  4u
  p16:gp0*1.1
  d16
#endif

  ; t1 evolution (with 1H and 15N CPD)
  4u pl8:f1 pl16:f3
  (p11 ph2):f1
  0.1u cpds1:f1 cpds3:f3

  (p3 ph3):f2
;  if "d0 > p14"
;  {
;  (center (d0) (p14:sp5 ph1):f2 )
;  } else {
  d0
;  }
  TAU do:f1 do:f3
  4u pl1:f1
  (center (p2 ph1):f1 (p14:sp3 ph5):f2 )
  TAU1
  p16:gp1*EA
  d16 pl2:f2
#ifdef ZZFILTER
  (p3 ph4):f2
  ; zz filter
  4u
  p16:gp3
  d16
  (p1 ph1)
#else
  (ralign (p1 ph1) (p3 ph4):f2 )
#endif /*ZZFILTER*/
  ; back-transfer
  4u
  p16:gp5
  d16
  DELTA4 pl0:f2
  (center (p2 ph1) (p8:sp13 ph1):f2 )
  4u
  p16:gp2
  DELTA1 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2 
  d1 do:f2 mc #0 to 2 
     F1EA(calgrad(EA), caldel(d0, +in0) & calph(ph3, +180) & calph(ph31, +180))
exit 
  

ph1=0
ph2=1
ph3=0 2
ph4=0 0 0 0 2 2 2 2
ph5=0 0 2 2
ph6=0
ph31=0 2 0 2 2 0 2 0


;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl6 : f1 channel - 3kHz water purge
;pl12: f2 channel - power level for CPD/BB decoupling
;sp3 : f2 channel - shaped pulse 180 degree (on resonance)
;sp5 : f2 channel - shaped pulse 180 degree (off resonance)
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p14: f2 channel - 180 degree shaped pulse
;p16: homospoil/gradient pulse
;p22: f3 channel - 180 degree high power pulse
;p28: f1 channel - trim pulse                        [1 msec]
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d20 : = d23
;d2: d2 = T : 26.6 or 53.2 msec
;     T (constant time period) = n/J(CC)
;cnst2: = J(XH)
;cnst21: CO chemical shift (offset, in ppm)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;in20: = in0
;nd0: 2
;ns: 4 * n
;ds: 32
;td1: number of experiments
;FnMODE: echo-antiecho
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    gp 1 : gp 2
;                         80 : 20.1    for C-13
;                         80 :  8.1    for N-15

;for z-only gradients:
;gpz1: 80%
;gpz2: 30.1% for C-13
;gpz3: 35% (zz filters)
;gpz4: 13% (180 refocusing)
;gpz5: 10% (180 refocusing)

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100




;$Id: hsqcctetgpsp,v 1.8.2.1 2015/03/03 11:21:23 ber Exp $
