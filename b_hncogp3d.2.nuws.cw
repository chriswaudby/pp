;b_hncogp3d.2.nuws.cw
;with NUWS in 15N dimension (for highly folded spectra)
;NB acquistion order 312 not 321
;Chris Waudby Jan 2018
;
;avance-version (15/03/12)
;best-HNCO
;3D sequence with
;   inverse correlation for triple resonance using multiple
;      inept transfer steps
;
;      F1(H) -> F3(N) -> F2(C=O,t1) -> F3(N,t2) -> F1(H,t3)
;
;on/off resonance Ca and C=O pulses using shaped pulse
;using shaped pulses for inversion and refocussing on f3
;phase sensitive (t1)
;phase sensitive using Echo/Antiecho (t2)
;using semi constant time in t2
;(use parameterset B_HNCOGP3D)
;
;P. Schanda, H. v. Melckebeke & B. Brutscher, 
;   J. Am. Chem. Soc. 128, 9042-9043 (2006)
;E. Lescop, P. Schanda & B. Brutscher, 
;   J. Magn. Reson.  187 163-169 (2007)
;(S. Grzesiek & A. Bax, J. Magn. Reson. 96, 432 - 440 (1992))
;(J. Schleucher, M. Sattler & C. Griesinger, 
;   Angew. Chem. Int. Ed. 32, 1489-1491 (1993))
;(L.E. Kay, G.Y. Xu & T. Yamazaki, J. Magn. Reson. A109, 
;   129-133 (1994))
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"d11=30m"

"d23=14.5m"
"d26=2.7m"

"p29=250u"

; NUWS prep
define loopcounter dsFlag
"dsFlag=1"

; number of complex points
"l3=td1/2"
"l6=td2/2"

#   ifdef CALC_SP
"p41=(bwfac25/(cnst55*cnst51*bf1))*1000000"
"spw25=plw1/((p41*90.0)/(p1*totrot25))*((p41*90.0)/(p1*totrot25))*(integfac25*integfac25)"
"spw27=plw1/((p41*90.0)/(p1*totrot27))*((p41*90.0)/(p1*totrot27))*(integfac27*integfac27)"
"spoal25=1"
"spoal27=0"

"p42=(bwfac26/(cnst55*cnst52*bf1))*1000000"
"spw26=plw1/((p42*90.0)/(p1*totrot26))*((p42*90.0)/(p1*totrot26))*(integfac26*integfac26)"
"spoal26=0.5"

"p43=(bwfac28/(cnst55*cnst53*bf1))*1000000"
"spw28=plw1/((p43*90.0)/(p1*totrot28))*((p43*90.0)/(p1*totrot28))*(integfac28*integfac28)"
"spw29=plw1/((p43*90.0)/(p1*totrot29))*((p43*90.0)/(p1*totrot29))*(integfac29*integfac29)"
"spoal28=1"
"spoal29=0"
#   endif /*CALC_SP*/


"d0=3u"
"d10=3u"
"d29=3u"
"d30=d23-p43-4u-p21*4/PI"

"in0=inf1/2"
"in10=inf2/2"

"FACTOR2=d30*10000000*2/td2"
"INCR2=FACTOR2/10000000"

"if ( INCR2 > in10 ) { in30 = in10; } else { in30 = INCR2; }"
"if ( INCR2 > in10 ) { in29 = 0; } else { in29=in10-INCR2; }"


"TAU=larger(p14,p44)"

"DELTA=d0*2+larger(TAU,p56)-p14"
"DELTA1=d26-p29-d16-p41*cnst41-larger(p42,p56)/2"
"DELTA2=d23-d26-p44-p16-d16-p14-d29"
"DELTA3=d26-p19-d16-p42/2"
"DELTA4=d26-p29-d16-p43*cnst43-larger(p42,p56)/2"
"DELTA5=p16+d16+de+8u" 
"DELTA6=d23-larger(p42,p57)/2"
"DELTA7=d23-larger(p42,p57)/2-p44-d26"
"DELTA8=d26-p14-d10"


"spoff2=0"
"spoff3=0"
"spoff5=bf2*(cnst22/1000000)-o2"
"spoff8=0"

"spoff25=bf1*(cnst54/1000000)-o1"
"spoff26=bf1*(cnst54/1000000)-o1"
"spoff27=bf1*(cnst54/1000000)-o1"
"spoff28=bf1*(cnst54/1000000)-o1"
"spoff29=bf1*(cnst54/1000000)-o1"
"spoff30=0"


aqseq 312


"acqt0=0"
baseopt_echo


1 d11 ze
  d11 pl26:f3 
2 d11 do:f3
3 d1
  50u UNBLKGRAD

  (p41:sp25 ph1)
  p29:gp3
  d16
  DELTA1
  (center (p42:sp26 ph1) (p56:sp39 ph1):f3 )
  DELTA1
  p29:gp3
  d16
  (p41:sp27 ph2):f1 

  p16:gp4
  d16 pl3:f3

  (p21 ph3):f3
  DELTA6
  (center (p14:sp3 ph1):f2 (p57:sp40 ph1):f3 )
  DELTA7
  (p44:sp30 ph1)
  d26 pl3:f3
  (p21 ph1):f3
  (p44:sp30 ph1)

  p16:gp5
  d16

  (p44:sp30 ph1)
  (p13:sp2 ph4):f2
  d0
  (center (p44:sp30 ph1) (p14:sp5 ph1):f2 (p56:sp39 ph7):f3 )
  d0
  4u
  (p14:sp3 ph1):f2
  DELTA
  (p14:sp5 ph1):f2
  4u
  (p13:sp8 ph1):f2

  p16:gp6
  d16 pl3:f3

  (p44:sp30 ph1)
  (p21 ph8):f3
  2u
  (p56:sp39 ph1):f3
  d10
  (p14:sp5 ph1):f2
  DELTA8
  (p44:sp30 ph1)
  DELTA2
  p16:gp1*EA
  d16 
  (p14:sp3 ph1):f2 
  d29
  (p56:sp39 ph7):f3
  d30
  2u pl3:f3
  (p43:sp28 ph1) 
  (p21 ph5):f3
  p19:gp7
  d16
  DELTA3
  (center (p42:sp26 ph1) (p57:sp40 ph1):f3 )
  DELTA3
  p19:gp7
  d16 pl3:f3
  (p21 ph6):f3

  (p43:sp29 ph2)
  p29:gp8
  d16
  DELTA4
  (center (p42:sp26 ph1) (p56:sp39 ph1):f3 )
  DELTA4
  p29:gp8
  d16
  (p43:sp28 ph1)
  DELTA5
  (p42:sp26 ph1)
  4u
  p16:gp2
  d16 pl26:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3 

; begin NUWS bit
  if "dsFlag==0" goto 10
  zd
  "dsFlag=0"
  goto 2  ; repeat following ds (without counting it as part of vclist)
10 4u

; repeat acquisition block according to schedule in vclist
  lo to 2 times c

; save data, reset scan counter
  4u do:f3
  d11 wr #0 if #0 zd 

; 13C looping (States-TPPI)
  1u ip4
  lo to 3 times 2
  1u id0
  1u ip4
  1u ip4
  lo to 3 times l3

; 15N looping (and NUWS incrementation)
  4u ivc
  1u rp4
  1u rd0
  1u igrad EA
  1u ip6
  1u ip6
  lo to 3 times 2
  1u id10
  1u id29
  1u dd30
  1u ip8
  1u ip8
  1u ip31
  1u ip31
  lo to 3 times l6



  d11 do:f3 mc #0 to 2 
     F1PH(calph(ph4, +90), caldel(d0, +in0)) 
     F2EA(calgrad(EA) & calph(ph6, +180), caldel(d10, +in10) & caldel(d29, +in29) & caldel(d30, -in30) & calph(ph8, +180) & calph(ph31, +180))
  TAU
exit


ph1=0
ph2=1 
ph3=0 0 0 0 2 2 2 2
ph4=0 2
ph5=0 0 2 2
ph6=1 1 3 3
ph7=0 0 0 0 2 2 2 2
ph8=0
ph31=0 2 2 0 2 0 0 2


;pl3 : f3 channel - power level for pulse (default)
;pl26: f3 channel - power level for CPD/BB low power decoupling
;sp2: f2 channel - shaped pulse  90 degree  (C=O on resonance)
;sp3: f2 channel - shaped pulse 180 degree  (C=O on resonance)
;sp5: f2 channel - shaped pulse 180 degree  (Ca off resonance)
;sp8: f2 channel - shaped pulse  90 degree  (C=O on resonance)
;                  for time reversed pulse
;sp25: f1 channel - shaped pulse  90 degree (Pc9_4_90.1000)
;sp26: f1 channel - shaped pulse 180 degree (Reburp.1000)
;sp27: f1 channel - shaped pulse  90 degree (Pc9_4_90.1000)
;                   for time reversed pulse
;sp28: f1 channel - shaped pulse  90 degree (Eburp2.1000)
;sp29: f1 channel - shaped pulse  90 degree (Eburp2tr.1000)
;                   for time reversed pulse
;sp30: f1 channel - shaped pulse 180 degree (Bip720,50,20.1)
;sp39: f3 channel - shaped pulse 180 degree (Bip720,50,20.1)
;sp40: f3 channel - shaped pulse 180 degree (Reburp.1000)
;p13: f2 channel -  90 degree shaped pulse
;p14: f2 channel - 180 degree shaped pulse
;p16: homospoil/gradient pulse                         [1 msec]
;p19: gradient pulse 2                                 [500 usec]
;p21: f3 channel -  90 degree high power pulse
;p29: gradient pulse 3                                 [250 usec]
;p41: f1 channel -  90 degree shaped pulse for excitation
;                      Pc9_4_90.1000             (2.2ms at 600.13 MHz)
;p42: f1 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000               (1.4ms at 600.13 MHz)
;p43: f1 channel -  90 degree shaped pulse for excitation
;                      Eburp2.1000/Eburp2tr.1000 (1.7ms at 600.13 MHz)
;p44: f1 channel - 180 degree shaped pulse for refocussing
;                      Bip720,50,20.1            (200us at 600.13 MHz)
;p56: f3 channel - 180 degree shaped pulse for inversion
;                      Bip720,50,20.1            (500us at 600.13 MHz)
;p57: f3 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000               (1.6ms at 600.13 MHz)
;d0 : incremented delay (F1 in 3D)                     [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d10: incremented delay (F2 in 3D)                     [3 usec]
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery
;d23: 1/(4J(NCO)                                       [14.5 msec]
;d26: 1/(4J(NH)                                        [2.7 msec]
;d29: incremented delay (F2 in 3D)                     [3 usec]
;d30: decremented delay (F2 in 3D) = d23-p43-4u-p21*4/PI
;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)
;cnst41: compensation of chemical shift evolution during p41
;           Pc9_4_90.1000: 0.529
;cnst43: compensation of chemical shift evolution during p43
;           Eburp2.1000: 0.5
;cnst51: scaling factor for p41 to compensate for transition region
;           Pc9_4_90.1000: 1.172
;cnst52: scaling factor for p42 to compensate for transition region
;           Reburp.1000: 1.426
;cnst53: scaling factor for p43 to compensate for transition region
;           Eburp2.1000: 1.000
;cnst54: H(N) chemical shift (offset, in ppm)
;cnst55: H(N) bandwidth (in ppm)
;o2p: CO chemical shift (cnst21)
;inf1: 1/SW(CO) = 2 * DW(CO)
;inf2: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(CO)) =  DW(CO)
;nd0: 2
;in10: 1/(2 * SW(N)) = DW(N)
;nd10: 2
;in29: = (1 - k2) * in10
;in30: = k2 * in10
;ns: 8 * n
;ds: >= 16
;aq: <= 50 msec
;td1: number of experiments in F1
;td2: number of experiments in F2
;FnMODE: States-TPPI (or TPPI) in F1
;FnMODE: echo-antiecho in F2
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 80%
;gpz2: 8.1%
;gpz3: 11%
;gpz4: 70%
;gpz5: 40%
;gpz6: 75%
;gpz7: 29%
;gpz8: 17%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.32
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam7: SMSQ10.50
;gpnam8: SMSQ10.32



                                          ;preprocessor-flags-start
;CALC_SP: for calculation of all bandselective Proton pulses based on cnst54 and cnst55
;             option -DCALC_SP (eda: ZGOPTNS)
                                          ;preprocessor-flags-end
										  


;$Id: b_hncogp3d.2,v 1.1.2.2 2015/03/12 17:07:07 ber Exp $
