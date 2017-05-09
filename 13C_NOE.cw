;13C_NOE
; 13C het-NOE measurement
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
"p0=p1*125/90" ; NOE saturation 125 degree pulse
"p4=p3*2"
"d3=1s/(cnst2*6.58)"
"d4=1s/(cnst2*4)"
"d11=30m"
"d12=20u"


"in0=inf2"
"d0=in0/2-p3*0.6366"

"DELTA1=d3-4u-p17-d16-larger(p1,p3)"
"DELTA2=d4-4u-p16-d16-larger(p1,p3)"

; loop counters for NOE during recycle period
"l0=1"
"d30=2.5m-p0*0.5"
"l4=d31/(p0+2*d30)"

"acqt0=0.6366*p1"
baseopt_echo


1 ze 
  d11 pl12:f2
2 d1 do:f2

  4u UNBLKGRAD
  4u pl1:f1 pl2:f2

  ; NOE buildup
  if "l0 %2 == 1"
  {
4 d30 
  (p0 ph1)
  d30
  lo to 4 times l4
  }
  else
  {
  d31
  }

  ; t1 evolution (with 1H CPD)
  d16 pl8:f1
  (p3 ph11):f2
  d0 cpd1:f1

  ; back-transfer
  4u do:f1
  p17:gp2
  d16
  DELTA1 pl1:f1
  (center (p2 ph1) (p4 ph1):f2 )
  4u
  DELTA1
  p17:gp2
  d16
  (p3 ph12):f2

  ; zz-filter
  4u
  p16:gp3
  d16

  ; second INEPT back-transfer
  (p1 ph1):f1
  4u
  p16:gp4
  d16
  DELTA2
  (center (p2 ph1) (p4 ph1):f2 )
  DELTA2
  p16:gp4
  d16 pl12:f2
  4u BLKGRAD

  ; acquisition
  go=2 ph31 cpd2:f2 
  d1 do:f2 mc #0 to 2
     F1QF(iu0) 
     F2PH(ip11, id0)
     F0()
exit 
  

ph1 =0
ph2 =1
ph11=0 2
ph12=1 1 3 3
ph31=0 2 2 0


;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
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
;p17: gradient pulse [300 usec]
;d0 : incremented delay (2D)
;d1 : relaxation delay; 1-5 * T1
;d31: saturation period
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(XH)
;cnst21: CO chemical shift (offset, in ppm)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(SW(X)) = DW(X)
;ns: 4 * n
;ds: 4
;td1: number of experiments
;FnMODE: states-tppi
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
