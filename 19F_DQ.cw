;19F DQ build-up (for measurement of CCR in A3 spin system)

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>



"d11=30m"
"d12=20u"

"p2=p1*2"

; for baseopt
"acqt0=-p1*2/3.1416"

1 ze 
2 d1 

  "DELTA=vd*0.5-p1"

  4u pl1:f1
  p1 ph1
  DELTA
  p2 ph1
  DELTA
  p1 ph1
  0.1u
  p1 ph2
  go=2 ph31 
  d1 mc #0 to 2 
     F1QF(ivd)

exit 
  

ph1 = 0 1 2 3
ph2 = {0 0 0 0}^1^2^3
ph31= {0 2 0 2}^3^2^1


;d1 : relaxation delay (excluding saturation time)
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;vdlist: relaxation time
;ns: 16 * n
;ds: 16


