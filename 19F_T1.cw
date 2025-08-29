;19F inversion recovery
;PJS updated for TS4.4 (onwards)

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>



"d11=30m"
"d12=20u"

"p2=p1*2"

; for baseopt
"acqt0=0"
baseopt_echo

define list<delay> t1_delay=<$VDLIST>	;UPDATED HERE

1 ze 
  
2 d1 

  "DELTA1=t1_delay"			;UPDATED HERE

  4u pl1:f1
  p2 ph1
  DELTA1
  p1 ph2
  go=2 ph31 
  d1 mc #0 to 2 
     F1QF(t1_delay.inc)			;UPDATED HERE

exit 
  

ph1 = 0 1 2 3
ph2 = 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph31= 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3


;d1 : relaxation delay (excluding saturation time)
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;vdlist: relaxation time
;ns: 4 * n
;ds: 4


