; zgfd.cw
; Calibration of water flipdown pulse
; Option for frequency-setting of shaped pulse


#include <Avance.incl>
#include <Delay.incl>

"d11=30m"
"d12=20u"

"acqt0=-p1*2/3.14159"

1 ze
2 d11
  d1 
  
# ifdef OFFRES
  10u fq=cnst20(bf ppm):f1
# endif /*OFFRES*/

  d12 pl0:f1
  (p29:sp11 ph2:r):f1
  10u pl1:f1

# ifdef OFFRES
  10u fq=0:f1
# endif /*OFFRES*/

  (p1 ph1)
  go=2 ph31
  d11 wr #0 
exit


ph1=0
ph2=2
ph31=0


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;sp11: f1 channel - shaped pulse  90 degree
;p1 : f1 channel -  90 degree high power pulse
;p29: f1 channel -  90 degree shaped pulse (flipback pulse)
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power-switching                        [20 usec]
;NS: 1 * n
;DS: 4

