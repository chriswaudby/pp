;19F broadband 1D T2 perfect-echo measurement
;with 1H decoupling

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<delay> t2delay = <$VDLIST>
/****************************/
/* Initialize loop counters */
/****************************/
"l1=0"


"d11=30m"
"d12=20u"

"p20=600u"
;"p21=12u"
;"spw20=plw1*pow(p1/p21,2)"

; for baseopt
;"acqt0=-p1*2/3.1416"
"acqt0=0"
"DELTA=1m"

1 ze 
  d11 pl12:f2
2 30m do:f2

  "DELTA=t2delay[l1]*0.25"

  ; purge
  ;20u pl11:f1
  ;(2mp ph11):f1
  ;20u
  ;(3mp ph12):f1

  ; d1
  d1

  ; 90
  (p20:sp20 ph1):f1  ; x (p2p)
  DELTA
  (p21:sp21 ph2):f1
  DELTA
  (p22:sp22 ph3):f1  ; y (90 unitary)
  DELTA
  (p21:sp21 ph4):f1
  DELTA

  go=2 ph31 cpd2:f2
  d11 do:f2 mc #0 to 2 
     F1QF(iu1)

exit 
 

ph1 =0 2
ph2 =0 0 1 1 2 2 3 3
ph3 =1
ph4 =0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1
     2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3
ph31=0 2 2 0 0 2 2 0 2 0 0 2 2 0 0 2

;pl12: f2 channel - power level for CPD/BB decoupling
;p16: homospoil/gradient pulse                       [0.5 msec]
;d1 : relaxation delay (excluding saturation time)
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence
;p20: 1000us
;spnam20: pulse_19F_Iz-Iy_0_15625_1000_100_0
;sp20: 16 kHz
;p21: 2000us
;spnam21: pulse_19F_180x_0_15625_2000_100_0
;sp21: 16 kHz
;p22: 1000us
;spnam22: pulse_19F_90x_0_15625_1000_100_0
;sp22: 16 kHz
;ns: 16 * n
;ds: 16


;for z-only gradients:
;gpz1: 41%

;use gradient files:   
;gpnam1: SMSQ10.100

