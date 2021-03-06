;hsqc_Nz_exsy.cw
;Nz exchange
;
;using vclist, 40 ms per repeat
;
;Chris Waudby, Nov 2017
;
;adapted from Farrow et al. (1994)
;
;run as pseudo-3D
;mixing times in seconds from vd-list

prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

;$DIM=3D
aqseq 312

define delay vdmin

"p2=p1*2"
"p22=p21*2"
"d4=1s/(cnst4*4)-p19-d16-larger(p21,p1)"
"d5=1s/(cnst4*4)-p19-d16-p11-4u-larger(p21,p1)"
"d11=30m"
"d12=20u"
"d15=50u"

"d0=0u"

"DELTA1=1s/(cnst4*4)-p19-d16-p21"

"in0=inf2"

"spoff5=bf1*(cnst19/1000000)-o1"
"d25=20m-0.5*p15"
"vdmin=8u+2*p19+2*d16"

"acqt0=0.6366*p1"

1 ze
  vdmin
  d11 pl16:f3
  d12  BLKGRAD
2 d1 do:f3
3 d12 pl3:f3

  ; calculate time for exchange delay
  10u
  ;"d21=vd-8u-2*p19-2*d16"
  10u

  ; Nz spoiler
  (p21 ph4):f3
  d15 UNBLKGRAD
  p19:gp1
  d16 pl1:f1

  ; first INEPT transfer
  (p1 ph4):f1
  p19:gp2
  d16
  d4
  (center (p2 ph4):f1 (p22 ph4):f3)
  d4
  p19:gp2
  d16

  (p1 ph1):f1
  4u pl0:f1
  (p11:sp1 ph4):f1
  4u pl1:f1

  p19:gp3
  d16

  (p21 ph5):f3

  p19:gp8
  d16

#ifdef LABEL_CN
  (center
    (p2 ph4):f1
    (p8:sp13 ph1):f2
    (DELTA1 d0 p22 ph6 DELTA1):f3
  )
#else
  (center
    ( p11:sp1 ph2 1u pl1 p2 ph4 1u p11:sp1 ph2 ):f1
    (DELTA1 d0 p22 ph6 DELTA1):f3
  )
#endif /*LABEL_CN*/

  p19:gp8
  d16 pl1:f1

  (p21 ph11):f3

  4u
  p19:gp4         ;Nz spoiler
  d16 BLKGRAD

;------------ Nz EXCHANGE PERIOD ------
77 d25
   (p15:sp5 ph4):f1
   d25
   lo to 77 times c

  4u UNBLKGRAD
  p19:gp5         ;Nz spoiler
  d16 pl1:f1

  (p21 ph4):f3

  p19:gp6
  d16
  d4
  (center (p2 ph4):f1 (p22 ph4):f3)
  d4
  p19:gp6
  d16

  (p21 ph1):f3
  4u pl0:f1
  (p11:sp1 ph7):f1
  4u pl1:f1
  (p1 ph8):f1

  p19:gp7
  d16
  d5 pl0:f1
  (p11:sp1 ph9):f1
  4u pl1:f1
  (center (p2 ph10):f1 (p22 ph4):f3 )
  4u pl0:f1
  (p11:sp1 ph9):f1
  d5 pl16:f3
  p19:gp7
  d16 BLKGRAD

  go=2 ph31 cpd3:f3
  d1 do:f3 mc #0 to 2
    F1QF(ivc)
    F2PH(ip5, id0)

exit

ph1= 1
ph2= 2
ph3= 3
ph4= 0
ph5= 0 2
ph6= 0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph11=1 1 3 3
ph7= 0
ph8= 0
ph9= 0 
ph10=2 
ph31=0 2 2 0 2 0 0 2 0 2 2 0 2 0 0 2 

;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl16 : f3 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p21 : f3 channel -  90 degree high power pulse
;p22 : f3 channel - 180 degree high power pulse
;p11 : f1 channel - water flipback [1 ms]
;cnst19 : 8.2 ppm
;spnam5 : Reburp.1000
;sp5 : f1 channel - power level for reburp
;p15 : f1 channel - reburp for CCR suppression during exchange period [1257 us (800 MHz)]
;p16 : watergate gradient
;p19 : spoiler gradient [500-900 us]
;sp1 : f1 channel - power level for water flipback
;spnam1 : sinc1.1000
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;d1 : relaxation delay; 1-5 * T1
;d16 : gradient recovery delay [200usec]
;cnst4: = J(NH)
;inf2: 1/SW(X) = 2 * DW(X)
;in0: 1/(1 * SW(X)) = 2 * DW(X)
;nd0: 1
;NS: 4 * n
;DS: 4 * n
;td1: number of experiments
;FnMODE: States-TPPI
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence

;gpz1: Eqm Nz crush [31 %]
;gpz2: 180 pair [8 %]
;gpz3: HzNz spoil [15 %]
;gpz4: Nz spoil [-20 %]
;gpz5: Nz spoil [-21 %]
;gpz6: 180 pair [16 %]
;gpz7: HzNz crush [51 %]
;gpz8: 180 pair [13.17 %]

;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam7: SMSQ10.100
;gpnam8: SMSQ10.100
