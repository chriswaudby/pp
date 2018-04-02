;adaptive_hsqcphpr.2d
; full 2D spectrum
;
;based on hsqcphpr
;avance-version (07/04/04)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive
;with decoupling during acquisition
;
;G. Bodenhausen & D.J. Ruben, Chem. Phys. Lett. 69, 185 (1980)
;
;$CLASS=HighRes
;$DIM=1D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Delay.incl>
#include <Grad.incl>


"p2=p1*2"
"p4=p3*2"

"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf1"
"d0=in0/2-2*0.6366*p3" ; initial evolution time for 90/180 phase correction
"d20=in0/2" ; actual first evolution time accounting for 13C pulse durations

"d4=1s/(cnst2*4)"
"DELTA=d4-larger(p2,p4)/2-p16-d16-4u"

"acqt0=0.6366*p1"

1 ze 
  d11 pl12:f2
2 d11 do:f2

  ; off-resonance presat
  30u pl9:f1
  30u fq=cnst21(bf hz):f1
  d1 cw:f1 ph1
  4u do:f1
  4u fq=0:f1
  d12 pl1:f1 pl2:f2 UNBLKGRAD

  ; purge 13C equilibrium magnetisation
  (p3 ph1):f2
  4u 
  p16:gp5
  d16
  4u

  ; INEPT
  p1 ph1
  4u
  p16:gp1
  d16
  DELTA 
  (center (p2 ph1) (p4 ph1):f2 )
  DELTA 
  p16:gp1
  d16
  4u
  (p1 ph2)

  ; zz-filter
  4u
  p16:gp2
  d16

  (p3 ph6):f2
  d0 ; 13C evolution
  (p3 ph7):f2

  ; zz-filter
  4u
  p16:gp3
  d16

  ; retro-INEPT
  (p1 ph2) 
  4u
  p16:gp4
  d16
  DELTA 
  (center (p2 ph1) (p4 ph1):f2 )
  DELTA pl12:f2 
  p16:gp4
  d16 BLKGRAD
  4u
  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2 F1PH(ip7, id0)
exit 
  

ph1=0
ph2=1
ph3=2
ph6=0 2
ph7=0 0 2 2
ph31=0 2 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p10 : f1 channel - 90 degree selective pulse [1000 usec]
;sp1 : f1 channel - 90 degree WFB (p10)

;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;cnst2: = J(XH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 1
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence

;for z-only gradients:
;gpz1: 7 %
;gpz2: 50 %
;gpz3: 35 %
;gpz4: 13 %
;gpz5: 17 %

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100

