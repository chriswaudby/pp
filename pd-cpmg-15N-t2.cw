;PD-CPMG (proton decoupled CPMG) for 15N T2 measurement
;Yuwen & Skrynnikov (2014) J Biomol NMR
;
;Chris Waudby, Dec 2018


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>



"p2=p1*2"
"p22=p21*2"

"d11=30m"
"d12=20u"

;"d4=1s/(cnst4*4)" ; 1/4J
"d4=2.25m" ; tau_a
"d5=2.75m" ; tau_b
"d21=d4-larger(p21,p22)-p32-d16"
"d22=d5-larger(p21,p22)-p33-d16"
"d23=d5-d0-p36-d16-larger(p1,p3)"
"d24=d5-p36-d16-larger(p1,p3)"
"d25=d4-larger(p21,p22)-p37-d16"
"d27=4u+de"
"acqt0=de"
baseopt_echo

"TAU=7.0218m-p22"
"TAU1=TAU*0.5"

define delay loopduration
"loopduration=4*7.0218m"

"in0=inf2"

aqseq 312


1 ze
  loopduration
  d11 pl16:f3
2 d11 do:f3

  ; purge water before d1
  4u UNBLKGRAD
  4u pl12:f1
  (6mp ph1)
  (3.7mp ph2)
  p39:gp9
  d16 pl1:f1
  (p1 ph1)
  p40:gp10
  d16
  4u BLKGRAD

  d1

  ; purge 15N magnetisation
  4u UNBLKGRAD
  4u pl3:f3
  (p21 ph1):f3
  p31:gp1
  d16

  ; begin main sequence - first INEPT
  (p1 ph1)
  p32:gp2
  d16
  d21
  (center (p2 ph1):f1 (p22 ph1):f3)
  d21
  p32:gp2
  d16
  (p1 ph2)

  ; zz filter
  p33:gp3
  d16

  ; second INEPT
  (p21 ph11):f3
  p33:gp3
  d16
  d22
  (center (p2 ph1):f1 (p22 ph1):f3)
  d22
  p33:gp3
  d16 pl13:f1 ; DIPSI power level

  ; turn on 1H decoupling *just* before relaxation block
  0.1u cpds1:f1

  ; T2 relaxation
77 TAU1
  (p22 ph1):f3
  TAU
  (p22 ph1):f3
  TAU2
  (p22 ph2):f3
  TAU
  (p22 ph4):f3
  TAU1
  lo to 77 times c

  ; deviation from published sequence - turn off decoupling during gradient purge
  (p21 ph3):f3
  4u do:f1
  p35:gp5
  d16

  ; restart decoupling for t1 evolution
  0.1u cpd1:f1
  
  ; t1 evolution
  (p21 ph12):f3
  d0
  
  ; first retro-INEPT
  4u do:f1
  p36:gp6*-1*EA
  d16
  d23 pl1:f1
  (center (p2 ph1):f1 (p22 ph13):f3 )
  d24
  p36:gp6*EA
  d16
  (ralign (p1 ph1):f1 (p21 ph14):f3 )

  ; second retro-INEPT
  p37:gp7
  d16
  d25
  (center (p2 ph1):f1 (p22 ph1):f3)
  d25
  p37:gp7
  d16
  (lalign (p1 ph2):f1 (p21 ph2):f3 )

  ; PEP
  p37:gp7
  d16
  d25
  (center (p2 ph1):f1 (p22 ph1):f3)
  d25
  p37:gp7
  d16
  (p1 ph4):f1

  ; spin-echo for gradient selection
  p38:gp8*-1
  d16
  d27
  (center (p2 ph1):f1 (p22 ph1):f3 )
  p38:gp8
  d16 pl16:f3
  4u BLKGRAD

  ; acquisition
  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 4
    F1QF(ivc)
    F2EA(igrad EA & ip14*2, id0 & ip12*2 & ip31*2)

  4u do:f3
exit
   

ph1=0 
ph2=1
ph3=2
ph4=3
ph11=0 2
ph12=1
ph13=0 0 1 1 2 2 3 3
ph14=0
ph31=0 2 2 0
  

;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl12: f1 channel - 12.4 kHz purge pulse
;pl13: f1 channel - 4.1 kHz DIPSI-2 decoupling (61 us)
;pl3 : f3 channel - power level for pulse (default)
;pl16: f3 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p31: gradient pulse [1 msec]
;p32: gradient pulse [500 usec]
;p33: gradient pulse [1 msec]
;p34: gradient pulse [500 usec]
;p35: gradient pulse [1 msec]
;p36: gradient pulse [1.25 msec]
;p37: gradient pulse [500 usec]
;p38: gradient pulse [125 usec]
;p39: gradient pulse [3.5 msec]
;p40: gradient pulse [2 msec]

;d1 : relaxation delay; 1-5 * T1
;d10 : incremented delay                             [3 usec]
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;vc : variable loop-coounter for T2 delay, taken from vc-list
;inf2: 1/SW(X) = 2 * DW(X)
;in10: 1/(2 * SW(X)) = DW(X)
;nd10: 2
;NS: 8 * n
;DS: >= 16
;td1: number of delays in vc-list
;FnMODE: QF in F1
;FnMODE: echo-antiecho in F2
;cpd1: DIPSI-2
;pcpd1: f1 channel - 4.1 kHz (61 us) DIPSI-2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 10%
;gpz2: 8%
;gpz3: 20%
;gpz4: 16%
;gpz5: 50%
;gpz6: 30%
;gpz7: 8%
;gpz8: 29.6%
;gpz9: 40%
;gpz10: 40%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SINE.20
;gpnam7: SMSQ10.100
;gpnam8: SINE.20
;gpnam9: SMSQ10.100
;gpnam10: SMSQ10.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: hsqct2etf3gpsi3d,v 1.5 2007/04/11 13:34:30 ber Exp $
