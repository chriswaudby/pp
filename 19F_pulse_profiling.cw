;19F Hahn-echo

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


define list<delay> t2delay = <$VDLIST>

"d11=30m"
"d12=20u"

; for baseopt
;"acqt0=-p1*2/3.1416"
"acqt0=0"
baseopt_echo

1 ze 
2 d11

  ; purge
  ;20u pl11:f1
  ;(2mp ph11):f1
  ;20u
  ;(3mp ph12):f1

  ; d1
  d1

  ; apply pulse
  4u pl20:f1
  (p20:sp20 ph1):f1

  t2delay*0.5

  (p21:sp21 ph2):f1

  t2delay*0.5

  go=2 ph31 
  d11 mc #0 to 2 
     F1QF(t2delay.inc)

exit 
  

ph1 =0 2 1 3
ph2 =0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph11=0
ph12=1
ph31=0 2 3 1 2 0 1 3


;p16: homospoil/gradient pulse                       [0.5 msec]
;d1 : relaxation delay (excluding saturation time)
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d18: saturation time
;pl8: f1 channel - power level for CEST saturation
;ns: 1 * n
;ds: 4


;for z-only gradients:
;gpz1: 41%

;use gradient files:   
;gpnam1: SMSQ10.100

