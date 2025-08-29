;19F CEST (scanning presat frequencies)

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


define list<frequency> fqlist = <$FQ1LIST>

"d11=30m"
"d12=20u"

; for baseopt
"acqt0=-p1*2/3.1416"

1 ze 
2 d1 

  ; CEST period
  4u pl8:f1 
  4u fq=fqlist:f1
  d18 cw:f1 ph11
  1u do:f1 

  ; purge
  4u UNBLKGRAD
  p16:gp1
  d16 pl1:f1 fq=0:f1
  4u BLKGRAD

/* ---------------------------------
; anti-ringing
; --------------------------------*/
 p1 ph1
 4u
 p1 ph2
 4u
 p1 ph3
;------------------------------------

  go=2 ph31 
  d1 mc #0 to 2 
     F1QF(calclist(fqlist,1))

exit 
  

ph1 =0
ph2 =2 0
ph3 =0 0 2 2 1 1 3 3
ph11=0 
ph31=0 2 2 0 1 3 3 1


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

