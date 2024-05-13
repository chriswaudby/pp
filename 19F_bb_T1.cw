;19F broadband 1D T1 measurement
;with 1H decoupling

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<delay> t1delay = <$VDLIST>
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

  "DELTA=t1delay[l1]-p16-d16-4u"

  ; purge
  ;20u pl11:f1
  ;(2mp ph11):f1
  ;20u
  ;(3mp ph12):f1

  ; d1
  d1

  4u UNBLKGRAD
  (p21:sp21 ph1):f1
  p16:gp1
  d16
  4u BLKGRAD
  DELTA
  ; 90 read-out
  (p20:sp20 ph1):f1

  go=2 ph31 cpd2:f2
  d11 do:f2 mc #0 to 2 
     F1QF(iu1)

exit 
 

ph1 =0 2 1 3
ph31=0 2 1 3

;pl12: f2 channel - power level for CPD/BB decoupling
;p16: homospoil/gradient pulse                       [0.5 msec]
;d1 : relaxation delay (excluding saturation time)
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence
;p20: 600us BURBOP_19F_90
;spnam20: BURBOP_19F_90
;sp20: 20 kHz
;ns: 1 * n
;ds: 4


;for z-only gradients:
;gpz1: 41%

;use gradient files:   
;gpnam1: SMSQ10.100

