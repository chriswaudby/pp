;imaging_zgesgp.cw
; 1D CSI
; td1 is hard-coded (look for size of r1d list)
;using baseopt
;avance-version (07/10/04)
;2D sequence
;water suppression using excitation sculpting with gradients
;T.-L. Hwang & A.J. Shaka, J. Magn. Reson.,
;   Series A 112 275-279 (1995)
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

; define list of gradient amplitudes, from -1 to 1, where final value is never reached
define list<grad_scalar, 384> r1d
;"td1=128"

"p2=p1*2"
"d12=20u"


"acqt0=de+4u+p19+d16+p1*0.6366"

1 ze
2 30m
#ifdef PRESAT
  d12 fq=cnst19(bf ppm):f1 ; off-resonance decoupling
  d12 pl9:f1 BLKGRAD
  d1 cw:f1
  4u do:f1
  4u pl1:f1
  d12 fq=0:f1 ; restore 1H frequency
#else
  d12 pl1:f1 BLKGRAD
  d1
  8u
#endif
  50u UNBLKGRAD
  p1 ph1

  ; imaging gradient
  p19:gp3*r1d
  d16
  4u BLKGRAD

  ; acquisition
  go=2 ph31
  30m mc #0 to 2 F1QF(r1d.inc)
  4u BLKGRAD
exit


ph1=0
ph31=0


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;p19: imaging gradient [200 usec]
;d1 : relaxation delay; 1-5 * T1
;d12: delay for power switching                             [20 usec]
;d16: delay for homospoil/gradient recovery
;NS: 1 * n, total number of scans: NS * TD0
;DS: 4

;gpz3: 17%

;use gradient files:
;gpnam3: SINE.32


