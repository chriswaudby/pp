;hsqcctetgpsppr.cw
;avance-version (21/09/15)
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
"d4=1s/(cnst2*4)"
"d11=30m"
"d12=20u"


"d0=3u"
"d20=d23-p16-d16-p14*1.5-4u-d12"

"in0=inf1/2"

"in20=in0"

"td1=tdmax(td1,d20*2,in20)"


"DELTA1=d4-larger(p2,p8)/2-p16-de-8u"
"DELTA2=d4-larger(p2,p8)/2-4u"
"DELTA3=d23-d0-p14/2-larger(p14,p22)-4u"
"DELTA4=d4-larger(p2,p8)/2-p1*2/PI-4u"


"spoffs3=0"
"spoffs5=bf2*(cnst21/1000000)-o2"
"spoffs13=0"


"acqt0=0"
baseopt_echo


1 ze 
  d11 pl12:f2
2 d11 do:f2
  4u BLKGRAD

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  30u fq=0:f1

3 (p1 ph1)
  DELTA2 pl0:f2
  4u
  (center (p2 ph1) (p8:sp13 ph6):f2 )
  4u
  DELTA2 pl2:f2 UNBLKGRAD

#   ifdef TRIMP
  (p28 ph1)
#   endif /* TRIMP */

  (p1 ph2) (p3 ph3):f2
  d0
  (center (p2 ph5) (p14:sp5 ph1):f2 (p22 ph1):f3 )
  4u
  DELTA3 pl0:f2
  (p14:sp3 ph4):f2
  d20
  p16:gp1*EA*-1
  d16 pl0:f2
  (p14:sp5 ph1):f2 
  4u
  d12 pl2:f2
  (ralign (p1 ph1) (p3 ph4):f2 )
  4u
  DELTA4 pl0:f2
  (center (p2 ph1) (p8:sp13 ph1):f2 )
  4u
  p16:gp2
  DELTA1 pl12:f2
  4u BLKGRAD
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 
     F1EA(calgrad(EA), caldel(d0, +in0) & caldel(d20, -in20) & calph(ph3, +180) & calph(ph6, +180) & calph(ph31, +180))
exit 
  

ph1=0
ph2=1
ph3=0 2
ph4=0 0 0 0 2 2 2 2
ph5=0 0 2 2
ph6=0
ph31=0 2 0 2 2 0 2 0
ph29=0

;pl0 : 0W
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
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
;d23: d23 = T : 13.3 or 26.6 msec
;     2T (constant time period) = n/J(CC)
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


;for z-only gradients:
;gpz1: 80%
;gpz2: 20.1% for C-13, 8.1% for N-15

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


                                          ;preprocessor-flags-start
;TRIMP: to use trimpulse p28@pl1 start experiment with
;          option -DTRIMP (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: $
