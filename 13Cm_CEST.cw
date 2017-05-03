;13C methyl CEST
;
;based on
;Bouvignies, G. & Kay, L. E.
;J. Biomol. NMR 53, 303â310 (2012).
;
;Adjusted delays for zero first-order phase correction
;
;sequence adapted from hsqcetgp.jk
;avance-version (07/04/04)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive using Echo/Antiecho-TPPI gradient selection
;with decoupling during acquisition
;using trim pulses in inept transfer
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<frequency> C13sat = <$FQ2LIST>


"p2=p1*2"
"p4=p3*2"
"d4=1s/(cnst2*4)"
"d5=1.216m" ; optimal tau_b transfer time
"d6=784u" ; optimal tau_c transfer time

"d11=30m"
"d12=20u"
"d13=4u"

"spoff1=cnst21-o1"

"in0=inf2/2"
"in10=inf2/2"
"l2=1" ; loop counter for CEST increments
"l3=td2/2"
"l6=td1"

"DELTA=d4-p16-d16-larger(p1,p3)-0.6366*p1"
"DELTA1=d5-p19-d16-larger(p1,p3)-0.6366*p3"
"d0=d6+in0/2-p19-d16-p3-0.6366*p3"
"d10=in0/2-p1-p3"
"DELTA2=d6-p19-d16-p1-0.6366*p3"
;"DELTA3=d4-p19-d16-p10-p1-4u-0.6366*p1"
;"DELTA4=d4-p19-d16-p10-p1-12u-de"
"DELTA3=d4-p19-d16-p1-4u-0.6366*p1"
"DELTA4=d4-p19-d16-p1-12u"
"acqt0=0"

aqseq 312


1 ze

  d11 pl1:f1 pl12:f2
2 d11 do:f2
  d11
  60u
3 90u
4 60u
5 10u

  ; temperature compensation following reference point
  10u
  if "l2==1"
  {
  10u
  10u pl8:f1
  10u LOCKH_ON
  d18 cpds1:f1
  10u do:f1
  10u LOCKH_OFF
  10u
  }
  10u pl1:f1

  ; purge any 1H magnetisation following acquisition
  50u UNBLKGRAD
  p16:gp0
  d16
  (p1 ph1)
  4u
  p16:gp0*0.71
  d16
  4u BLKGRAD

  ; off-resonance presat
  30u pl9:f1
  30u fq=cnst21(bf hz):f1
  d1 cw:f1 ph1
  30u do:f1
  30u pl1:f1
  30u fq=0:f1

  ; purge equilibrium 13C
  30u UNBLKGRAD
  4u pl2:f2
  (p3 ph1):f2
  p16:gp0
  d16

  ; begin main sequence
  (p1 ph1)
  p16:gp1
  d16
  DELTA
  (center (p2 ph1) (p4 ph1):f2 )
  DELTA
  p16:gp1
  d16
  (p1 ph2)

  ; zz purge
  p16:gp2
  d16

  (p3 ph11):f2
  p19:gp3
  d16
  DELTA1
  (center (p2 ph1) (p4 ph1):f2 )
  DELTA1
  p19:gp3
  d16
  (p3 ph2):f2

  (p1 ph1)
  4u
  p16:gp4
  d16

  ; CEST period
  10u
  if "l2==1" goto 77
  10u

  10u pl8:f1 pl18:f2
  10u LOCKH_ON
  10u C13sat:f2
  d18 cpds1:f1 cw:f2 ph1
  10u do:f1 do:f2
  10u LOCKH_OFF
  10u fq=0:f2

  10u
  77 10u
  10u

  5u pl1:f1 pl2:f2

  5u
  p16:gp5  ; cleaning gradient
  d16

  ; t1 evolution
  (p3 ph12):f2
  p19:gp6
  d16
  d0
  (p2 ph1)
  d10
  (p4 ph1):f2
  DELTA2
  p19:gp6
  d16
  (p3 ph2):f2

  ; zz purge
  p16:gp7
  d16

  ; final inept
  (p1 ph1)
  p19:gp8
  d16
  DELTA3
  ;(p10:sp1 ph3):f1
  4u
  (center (p2 ph1) (p4 ph1):f2 )
  4u
  ;(p10:sp1 ph3):f1
  DELTA4
  p19:gp8
  d16
  4u BLKGRAD
  4u pl12:f2

  go=2 ph31 cpd2:f2
  d11 do:f2
  d11 wr #0 if #0 zd

  ; loop over CEST frequencies
  30u C13sat.inc
  30u iu2
  lo to 3 times l6

  ; Echo/Anti-echo
  30u ip12
  30u ru2
  30u C13sat.res
  lo to 4 times 2

  ; t1 evolution
  30u id0
  30u id10
  lo to 5 times l3

  ;1m do:f2
  1m BLKGRAD




exit


ph1=0
ph2=1
ph3=2
ph4=3
ph11=0 2
ph12=0 0 2 2
ph31=0 2 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p10 : f1 channel - 90 degree selective pulse [1000 usec]
;sp1 : f1 channel - 90 degree WFB (p10)
;p16: homospoil/gradient pulse [1000 usec]
;p19: second homospoil/gradient pulse [300 usec]
;p22: f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                  [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(XH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 1 * n
;DS: >= 16
;td1: number of experiments
;FnMODE: states-tppi

;pcpd1: 1H 90 for CPD (4.5 kHz at 900 MHz = 55.5 us)
;pl8: 1H decoupling power (4.5 kHz at 900 MHz = 55.5 us)

;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence
;CPDPRG2: waltz16

;d18: relaxation time
;pl18: 13C CEST saturation power
;l6: number of CEST frequencies (including reference) (td1)
;l3: number of complex points (td2 / 2)

;gradients
;p16: 1000u
;p19: 300u
;
;gpz0: 15%
;gpz1: 10%
;gpz2: 24%
;gpz3: -16% [p19]
;gpz4: -32%
;gpz5: 14%
;gpz6: -24% [p19]
;gpz7: 5%
;gpz8: 30% [p19]

;use gradient files:
;gpnam0: SMSQ10.100
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SINE.10
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SINE.10
;gpnam7: SMSQ10.100
;gpnam8: SINE.10


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: hsqcetgp,v 1.4 2007/04/11 13:34:30 ber Exp $
