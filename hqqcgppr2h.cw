; simultaneous HQQC and HDQQC 
; with 2h decoupling during t1
; set td1 = 21 * actual td1
; requires post-processing for receiver phase
; Chris Waudby Jul 2020
;

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf1"
"d0=in0/2"

"DELTA1=d2-0.6366*p1-0.5*larger(p2,p3)"
"DELTA2=d2-0.5*larger(p2,p3)-0.5*larger(p1,p4)-p30-8u"
"DELTA3=d2-0.5*p1-0.5*larger(p2,p3)-p30-8u"
"DELTA4=d2-0.5*larger(p2,p3)-d12-4u"

"acqt0=0"
baseopt_echo

"l1=0"
"l3=td1/21"

1 ze 
  d11 LOCKDEC_ON
  50u LOCKH_ON
  d11 H2_PULSE

  d11 pl12:f2
  d11 pl17:f4
2 d11 do:f2
  4u BLKGRAD
  d11 H2_LOCK ; lock-in over d1
  6m LOCKH_OFF

  ; relaxation delay
# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  30u fq=0:f1

  ; turn lock hold on for 2H pulsing during main sequence
  d11 LOCKH_ON
  d11 H2_PULSE
  50u UNBLKGRAD
  d12 pl1:f1 pl2:f2 pl17:f4

  ; purge Cz
  (p3 ph1):f2
  4u
  p16:gp1
  d16*2 

  ; begin main sequence
  (p1 ph11):f1
  DELTA1 ; 1/2J
  (center (p2 ph11):f1 (p3 ph12):f2 )
  DELTA2

  ;switch on 2H decoupling just before t1
  (p30 ph21):f4
  4u
  4u cpds4:f4 ph20

  ; t1 evolution
  (center (p1 ph11):f1 (p4 ph1):f2 )
  d0
  p1 ph1

  ; 2H decoupling off and flip back
  4u do:f4
  4u
  (p30 ph23):f4

  ; back-transfer
  DELTA3
  (center (p2 ph1):f1 (p3 ph13):f2 )
  DELTA4
  4u BLKGRAD
  d12 pl12:f2

  go=2 ph31 cpd2:f2
  d12 do:f2
  d11 wr #0 if #0 zd
  4u iu1
  4u ip12*2
  if "l1 % 3 == 0"
  {
    4u ip11
  }
  lo to 2 times 21
  4u ru1
  4u id0
  4u ip13*2
  lo to 2 times l3

  d11 H2_LOCK
  d11 LOCKH_OFF
  d11 LOCKDEC_OFF
exit 
  
  
ph1=0 
ph2=1
ph3=2
ph4=3
ph11=(7) 0
ph12=(6) 0 3 
ph13=0
ph20=0
ph21=1
ph23=3
ph29=0
ph31=0 2


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;pl17: f4 channel - power level for 2H CPD decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p16: homospoil/gradient pulse
;p22 : f3 channel - 180 degree high power pulse
;p30: f4 channel - 90 degree pulse at pl17 (CPD 90)
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 21 * n
;DS: 21
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz2: 80%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100

