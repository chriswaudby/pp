;B1 calibration
;for 1H CEST with saturation during zz period
;based on Sekhar et al. PNAS 113, E2794-801 (2016)
;
;Chris Waudby (c.waudby@ucl.ac.uk)


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>



"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"
"d13=4u"
"d26=1s/(cnst4*4)"


"d18=1u"
;"in18=inf1" ; nutation time increment
;"l3=td1"


"DELTA=d26-p21-p16-d16-4u"
"DELTA1=d26-p21-p16-d16-8u-p11"
"DELTA2=d26-p21-p16-d16-12u-p11"
"acqt0=0" ; baseopt

; initial delay for half dwell phase correction
"d0=in0/2-1.27324*p21-p2-12u"




1 ze
  d11 
2 d1 
  d12 pl1:f1 pl3:f3

  50u UNBLKGRAD 	; crush equilibrium Nz (p16, gp4)
  (p21 ph10):f3
  p16:gp1
  d16*2

  (p10:sp10 ph13:r):f1   ; H2O flip-down (sp10)
  4u pl1:f1
  (p1 ph10):f1		; start INEPT
  4u
  p16:gp2
  d16
  DELTA
  (center (p2 ph11):f1 (p22 ph10):f3 )
  DELTA
  4u
  p16:gp2
  d16
  (p1 ph11):f1

  ; zz purge
  4u
  p16:gp3
  d16

  ; 15N evolution

  (p21 ph1):f3
  3u
;  d0*0.5 gron4
;  3u groff
  (p2 ph11):f1
  3u
  ;d0*0.5 gron4*-1
  ;3u groff
  (p21 ph2):f3

;goto 999  ;optimise sp10

  ; zz purge
  4u
  p16:gp5
  d16

  ; CEST period
  4u pl8:f1
  4u fq=cnst19(bf ppm):f1

  4u LOCKH_ON
  d18 cw:f1 ph10
  4u do:f1
  4u LOCKH_OFF

  4u pl1:f1
  4u fq=0:f1

  ; zz purge
  4u
  p16:gp6
  d16
  4u pl0:f1

  ; retro-INEPT
  (p11:sp11 ph12):f1   ; H2O flip-down (sp11)
  4u pl1:f1
  (p1 ph10)
;goto 999 ; optimise sp11
  4u
  p16:gp7
  d16
  DELTA1
  (p11:sp12 ph12):f1
  4u pl1:f1
  (center (p2 ph10):f1 (p22 ph10):f3)
  4u pl0:f1
  (p11:sp12 ph12):f1
  DELTA2
  p16:gp7
  d16
999  4u BLKGRAD
  4u 

  ; acquisition (without decoupling)
  go=2 ph31 

  30m mc #0 to 2 F0(zd)


exit


ph1=  0 2
ph2=  0
ph10= 0
ph11= 1
ph12= 2
ph13= 3
ph31= 0 2

;pl0 : 120 dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl8 : f1 channel - power level for CEST saturation
;pl16: f3 channel - power level for CPD/BB decoupling
;sp10: 90 deg Squa1.1000 water flip-down pulse
;sp11: 90 deg Squa1.1000 water flip-down pulse
;sp12: 90 deg Squa1.1000 water flip-down pulse
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p10: 1000u 90 deg soft rectangular water flip-back pulse
;p11: 1000u 90 deg soft rectangular water flip-back pulse
;p16: homospoil/gradient pulse [1 ms]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;d18: saturation time (incremented)
;d26 : 1/(4J(YH))
;cnst4: = J(YH)
;cnst19: saturation frequency for B1 calibration
;inf2: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 2 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI (or TPPI)
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 37%
;gpz2: 23%
;gpz3: 49%
;gpz4: 0.5%
;gpz5: 31%
;gpz6: 57%
;gpz7: 11%

;use gradient files:
;gpnam1: SQUA.100
;gpnam2: SQUA.100
;gpnam3: SQUA.100
;gpnam5: SQUA.100
;gpnam6: SQUA.100
;gpnam7: SQUA.100

                                          ;preprocessor-flags-start
                                          ;preprocessor-flags-end
