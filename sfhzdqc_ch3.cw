;methyl-SOFAST-H(Z/D)QC
;with multiplet filter
;run as 2D with td = 6 * actual td
;aqmod QF - sw WILL be correctly calculated
;
;Chris Waudby, 9/10/18
;
;Added option for off-resonance presat (e.g. to suppress urea signal), 21/6/15


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"d11=30m"
"d12=20u"
"d13=4u"
"d21=1s/(cnst2*2)"
"d22=1s/(cnst2*8)" 
"p4=p3*2"

"in0=inf1"
"d0=in0/2"


"DELTA1=d21-p16-d16*2-p39*cnst39"
"DELTA2=d21-p16-d16*2-4u-de"
"DELTA3=d22+d16-p40*0.5"
"DELTA4=d22-p3*1.6366"
"acqt0=de"


# ifdef OFFRES_PRESAT
  "TAU=d1-10m-60u-d12*2-d13"
# else
  "TAU=d1-10m"
# endif /*OFFRES_PRESAT*/


"spoff23=bf1*(cnst19/1000000)-o1"
"spoff24=bf1*(cnst19/1000000)-o1"


"l0=0"


1 ze 
  d11 pl12:f2
2 10m do:f2

# ifdef OFFRES_PRESAT
  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  TAU cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  30u fq=0:f1
# else
  TAU
# endif /*OFFRES_PRESAT*/

3 d12 pl2:f2
  50u UNBLKGRAD

  ; purge Cz
  (p3 ph1):f2
  p16:gp2
  d16

  ; begin main sequence
  (p39:sp23 ph13):f1
  p16:gp1
  d16
  DELTA1

  if "l0 % 4 == 0"
     {
  (lalign (DELTA3 p40:sp24 ph1) (d16 p3 ph11 DELTA4 p4 ph1 d0 DELTA4 p3 ph12 d16):f2)
     }
  if "l0 % 4 == 1"
     {
  (lalign (DELTA3 p40:sp24 ph1) (d16 p3 ph21 DELTA4 d0 p4 ph1 DELTA4 p3 ph22 d16):f2)
     }
  if "l0 % 4 == 2"
     {
  (lalign (DELTA3 d0 p40:sp24 ph1) (d16 p3 ph11 DELTA4 d0 p4 ph1 DELTA4 p3 ph12 d16):f2)
     }
  if "l0 % 4 == 3"
     {
  (lalign (DELTA3 d0 p40:sp24 ph1) (d16 p3 ph21 DELTA4 p4 ph1 d0 DELTA4 p3 ph22 d16):f2)
     }

  DELTA2
  p16:gp1
  d16 pl12:f2
  4u BLKGRAD

  go=2 ph31 cpd2:f2
  10m do:f2 mc #0 to 2 
     F1I(ip11*2 & ip21*2, 3, iu0, 4)
     F1QF(id0)

exit 
  

ph1=0 
ph11=(6) 0 ;2 4
ph21=(6) 3 ;5 1
ph12=0 2
ph22=2 0
ph13=0 0 2 2
ph29=0
ph31=0 2 2 0


;pl2 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling (low power)
;sp23: f1 channel - shaped pulse 120 degree (Pc9_4_120.1000)
;sp24: f1 channel - shaped pulse 180 degree (Reburp.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p3:  f2 channel -  90 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000 
;d0 : incremented delay (2D) 
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(CH) [125 Hz for methyls]
;cnst19: H(met) chemical shift (offset, in ppm)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_120.1000: 0.529
;inf1: 1/SW(C) = 2 * DW(C)
;in0: 1/ SW(C) = 2 * DW(C)
;nd0: 1
;NS: 2 * n
;DS: 16
;aq: <= 50 msec
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEC
;cpd2: decoupling according to sequence defined by cpdprg2: garp4
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:	gp 1 : gp 2
;			  11 :    7


;for z-only gradients:
;gpz1: 11%
;gpz2:  7%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100




;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1


