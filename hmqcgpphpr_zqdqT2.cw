;methyl ZQ and DQ T2 measurement using HMQC with multiplet filter
;Chris Waudby, July 2020
;
;set td2 = 6 * desired td (3 step cycle for ZQ/DQ selection + 2 step cycle for multiplet suppression)
;
;avance-version (07/04/04)
;HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive
;with decoupling during acquisition
;
;A. Bax, R.H. Griffey & B.L. Hawkins, J. Magn. Reson. 55, 301 (1983)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d3=1s/(cnst2*8)"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf2"
"d0=2*0.63662*p3 + in0/2"

; loop counter for ZQ/DQ and E/AE blocks
"l1 = 0"

"acqt0=0"
baseopt_echo

aqseq 321

1 ze 
  d11 pl12:f2
2 d11 do:f2
  ; purge before d1
  20u pl6:f1
  (2mp ph1):f1
  (3mp ph2):f1

  4u BLKGRAD

  ; relaxation period, with off-res presat
  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1 pl2:f2
  30u fq=0:f1
  50u UNBLKGRAD

  (p3 ph1):f2    ; crush eq'm 13C magnetisation
  d13
  p16:gp1
  d16

  ; start main sequence
  (p1 ph1):f1  ; INEPT
  "DELTA1 = d2 - p1*0.6366"
  DELTA1

  if "l1 % 4 == 0"
    {
  (p3 ph11):f2
  "DELTA = d3 - p3"
  DELTA
  (p3*2 ph1):f2
  "DELTA = vd*0.5 + d0*0.5 - p3 - p1*2 - p19 - d16"
  DELTA
  p19:gp2
  d16
  (p1 ph2):f1
  (p2 ph1):f1
  (p1 ph2):f1
  "DELTA = d3 + d0*0.5 - p3 - p1*2"
  DELTA
  (p3*2 ph1):f2
  p19:gp2
  d16
  "DELTA = vd*0.5 - p3 - p19 - d16"
  DELTA
  (p3 ph1):f2
    }

  if "l1 % 4 == 1"
    {
  (p3 ph11):f2
  (p3*2 ph1):f2
  "DELTA = d3 - p3 + vd*0.5 + d0*0.5 - p3 - p1*2 - p19 - d16"
  DELTA
  p19:gp2
  d16
  (p1 ph2):f1
  (p2 ph1):f1
  (p1 ph2):f1
  "DELTA = d0*0.5 - p1*2"
  DELTA
  (p3*2 ph1):f2
  p19:gp2
  d16
  "DELTA = d3 - p3 + vd*0.5 - p3 - p19 - d16"
  DELTA
  (p3 ph1):f2
    }

  if "l1 % 4 == 2"
    {
  (p3 ph11):f2
  "DELTA = vd*0.5 - p3 - p19 - d16"
  DELTA
  p19:gp2
  d16
  (p3*2 ph1):f2
  "DELTA = d3 + d0*0.5 - p3 - p1*2"
  DELTA
  (p1 ph2):f1
  (p2 ph1):f1
  (p1 ph2):f1
  p19:gp2
  d16
  "DELTA = vd*0.5 + d0*0.5 - p3 - p1*2 - p19 - d16"
  DELTA
  (p3*2 ph1):f2
  "DELTA = d3 - p3"
  DELTA
  (p3 ph1):f2
    }

   if "l1 % 4 == 3"
    {
  (p3 ph11):f2
  "DELTA = vd*0.5 - p3 - p19 - d16 + d3 - p3"
  DELTA
  p19:gp2
  d16
  (p3*2 ph1):f2
  "DELTA = d0*0.5 - p1*2"
  DELTA
  (p1 ph2):f1
  (p2 ph1):f1
  (p1 ph2):f1
  p19:gp2
  d16
  "DELTA = vd*0.5 + d0*0.5 - p3 - p1*2 - p19 - d16 + d3 - p3"
  DELTA
  (p3*2 ph1):f2
  (p3 ph1):f2
    }

   ; back-transfer
  "DELTA2 = d2 - d12 - 4u"
  DELTA2
  d12 pl12:f2
  4u BLKGRAD

  ; acquisition
  go=2 ph31 cpd2:f2 

  d11 do:f2 mc #0 to 2
    F2I(ip11, 3, iu1, 2)
    F2EA(iu1, id0)
    F1QF(rd0 & rp11 & ru1 & ivd)

  4u BLKGRAD
exit 
  

ph1= 0 
ph2= 1 
ph11= (3) 0
ph29=0
ph31=0


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse
;p19: gradient pulse [300 usec]
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d3 : 1/(8J)CH
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;l0: number of repeats for entire experiment
;NS: 1 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz2: 7%
;gpz3: -40%
;gpz4: 29%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.32
;gpnam4: SMSQ10.32
