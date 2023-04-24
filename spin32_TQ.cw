;TQ filtered relaxation for measurement of spin 3/2 systems
;pseudo2D experiment using vdlist
;
;https://www.sciencedirect.com/science/article/abs/pii/002223649090204M
;https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3634872/#FD1
;https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2896046/

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
  p1 ph2
  4u
  p1 ph3
  go=2 ph31 
  d1 mc #0 to 2 
     F1QF(ivd)

exit 
  

ph1 = (12) 1 3 5 7 9 11
ph2 = (12) 4 6 8 10 12 2
ph3 = 1
ph31= 0 2


;d1 : relaxation delay                               [1 sec]
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;vdlist: relaxation time
;ns: 6 * n
;ds: 6 * n


