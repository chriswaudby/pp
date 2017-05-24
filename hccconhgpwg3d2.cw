;With low-power presaturation to alleviate radiation damping
;
;hccconhgpwg3d2
;avance-version (08/11/06)
;HCCCONH
;3D sequence with
;   inverse correlation for triple resonance using 
;      dipsi2 and inept transfer steps 
;
;      F1(H,t1) -> F2(Caliph. -> Ca) -> F2(C=O) -> F3(N,t2) -> F1(H,t3)
;
;on/off resonance Ca and C=O pulses using shaped pulse
;phase sensitive (t1)
;phase sensitive (t2)
;using semi-constant time in t1
;using constant time in t2
;water suppression using watergate sequence
;(use parameterset HCCCONHGPWG3D2)
;
;G.T. Montelione, B.A. Lyons, S.D. Emerson & M. Tashiro, 
;   J. Am. Chem. Soc. 114, 10974-75 (1992)
;S. Grzesiek, J. Anglister & A. Bax, J. Magn. Reson. 101 B, 114-9 (1993)
;B.A. Lyons & G.T. Montelione, J. Magn. Reson. 101 B, 206-9 (1993)
;T.M. Logan, E.T. Olejniczak, R.X. Xu & S.W. Fesik, 
;   J. Biomol. NMR 3, 225-31 (1993)
;R.T. Clowes, W. Boucher, C.H. Hardman, P.J. Domaille & E.D. Laue,
;   J. Biomol. NMR 3, 349-354 (1993)
;T. Carlomagno, M. Maurer, M. Sattler, M.G. Schwendinger, S.J. Glaser
;   & C. Griesinger, J. Biomol. NMR 8, 161-170 (1996)
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include<Avance.incl>
#include<Grad.incl>
#include<Delay.incl>


"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"
"d13=4u"

"d3=1.1m"
"d4=1.7m"
"d21=3.6m"
"d22=4.4m"
"d23=12.4m"
"d25=5.5m"
"d26=2.3m"


"d0=d4"
"d10=d23/2-p14/2"
"d20=3u"
"d28=d4+p14+d20"
"d29=d23/2-p14/2-p26-d25-4u"
"d30=d23/2-p14/2"

"in0=inf1/2"
"in10=inf2/4"

"FACTOR1=trunc(d28*10000000*2/(td1-2))"
"in28=FACTOR1/10000000"

"in20=in0-in28"
"in29=in10"
"in30=in10"


"l1=(d15/(p9*115.112))+0.5"


"DELTA1=d23+4u-d22"
"DELTA2=d26-p16-d16-p11-12u"


"spoff2=0"
"spoff3=0"
"spoff5=bf2*((cnst21-cnst23)/1000000)"
"spoff7=bf2*((cnst22-cnst21)/1000000)"
"spoff8=0"
"spoff9=bf2*((cnst22-cnst23)/1000000)"


aqseq 321


1 ze

  "if ( in28 > in0 ) { in28=in0; in20= 0; }"

  d11 pl16:f3
2 d11 do:f3
3 d11 fq=cnst23(bf ppm):f2
  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  50u UNBLKGRAD
  d12 pl1:f1 pl0:f2 pl3:f3

  (p1 ph3):f1
  d0
  (p14:sp3 ph1):f2
  d20
  (p2 ph1):f1
  d28
  (p1 ph2):f1

  4u
  p16:gp1
  d16

  (p13:sp2 ph1):f2
  d3
  (center (p2 ph1):f1 (p14:sp3 ph1):f2 )
  d3
  (p13:sp8 ph2):f2

  4u
  d12 pl15:f2

						;begin DIPSI2
7 (p9*3.556 ph23):f2
  (p9*4.556 ph25):f2
  (p9*3.222 ph23):f2
  (p9*3.167 ph25):f2
  (p9*0.333 ph23):f2
  (p9*2.722 ph25):f2
  (p9*4.167 ph23):f2
  (p9*2.944 ph25):f2
  (p9*4.111 ph23):f2
  
  (p9*3.556 ph25):f2
  (p9*4.556 ph23):f2
  (p9*3.222 ph25):f2
  (p9*3.167 ph23):f2
  (p9*0.333 ph25):f2
  (p9*2.722 ph23):f2
  (p9*4.167 ph25):f2
  (p9*2.944 ph23):f2
  (p9*4.111 ph25):f2

  (p9*3.556 ph25):f2
  (p9*4.556 ph23):f2
  (p9*3.222 ph25):f2
  (p9*3.167 ph23):f2
  (p9*0.333 ph25):f2
  (p9*2.722 ph23):f2
  (p9*4.167 ph25):f2
  (p9*2.944 ph23):f2
  (p9*4.111 ph25):f2

  (p9*3.556 ph23):f2
  (p9*4.556 ph25):f2
  (p9*3.222 ph23):f2
  (p9*3.167 ph25):f2
  (p9*0.333 ph23):f2
  (p9*2.722 ph25):f2
  (p9*4.167 ph23):f2
  (p9*2.944 ph25):f2
  (p9*4.111 ph23):f2
  lo to 7 times l1
						;end DIPSI2

  d12 pl19:f1 pl0:f2
  (p26 ph1):f1
  d12 cpds1:f1 ph2

  (p13:sp2 ph4):f2
  4u
  (p14:sp5 ph1):f2
  d21 
  (p24:sp9 ph1):f2
  4u
  (p14:sp5 ph1):f2
  d21
  (p13:sp8 ph2):f2

  4u
  d12 fq=cnst21(bf ppm):f2

  (p13:sp2 ph5):f2
  d22
  (p14:sp7 ph1):f2
  DELTA1 pl3:f3
  (center (p14:sp3 ph1):f2 (p22 ph8):f3 )
  d23
  (p14:sp7 ph1):f2
  4u
  (p13:sp8 ph1):f2

  4u
  4u do:f1
  (p26 ph10):f1
  4u
  p16:gp2
  d16
  (p26 ph2):f1
  4u cpds1:f1 ph1

  (p21 ph6):f3
  d30
  (p14:sp7 ph1):f2
  d30
  (center (p14:sp3 ph1):f2 (p22 ph1):f3 )
  d10
  (p14:sp7 ph1):f2
  d29
  4u do:f1
  (p26 ph9):f1
  d25
  (p21 ph1):f3

  p16:gp3
  d16 
  (p11:sp1 ph7):f1
  4u
  4u pl1:f1

  (p1 ph1) 
  4u
  p16:gp4
  d16
  DELTA2 
  (p11:sp1 ph7):f1
  4u
  4u pl1:f1
  (center (p2 ph1) (p22 ph1):f3 )
  4u
  (p11:sp1 ph7):f1
  4u
  DELTA2
  p16:gp4
  d16 pl16:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
     F1PH(rd10 & rd29 & rd30 & rp6 & ip3, id0 & id20 & dd28 & dp3)
     F2PH(dp6, id10 & id29 & dd30)
exit


ph1=0
ph2=1
ph3=0 0 0 0 2 2 2 2
ph4=0 0 2 2
ph5=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph6=0 2
ph7=2
ph8=0
ph9=3
ph10=2
ph23=0
ph25=2
ph29=0
ph31=0 2 2 0 2 0 0 2 2 0 0 2 0 2 2 0


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl15: f2 channel - power level for TOCSY-spinlock
;pl16: f3 channel - power level for CPD/BB decoupling
;pl19: f1 channel - power level for CPD/BB decoupling
;sp1: f1 channel - shaped pulse  90 degree  (H2O on resonance)
;sp2: f2 channel - shaped pulse  90 degree  (C=O on resonance)
;sp3: f2 channel - shaped pulse 180 degree  (on resonance)
;sp5: f2 channel - shaped pulse 180 degree  (C=O off resonance)
;sp7: f2 channel - shaped pulse 180 degree  (Ca off resonance)
;sp8: f2 channel - shaped pulse  90 degree  (on resonance)
;                  for time reversed pulse
;sp9: f2 channel - shaped pulse 180 degree  (Ca on resonance)
;     sp9 might require higher selectivity than sp3
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p9 : f2 channel -  90 degree low power pulse
;p11: f1 channel -  90 degree shaped pulse             [1 msec]
;p13: f2 channel -  90 degree shaped pulse
;p14: f2 channel - 180 degree shaped pulse
;p16: homospoil/gradient pulse                         [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p24: f2 channel - 180 degree shaped pulse (sp9)
;p26: f1 channel -  90 degree pulse at pl19
;d0 : incremented delay (F1 in 3D) = d4
;d1 : relaxation delay; 1-5 * T1
;d3 : 1/(6J(CH)                                        [1.1 msec]
;d4:  1/(4J(CH)                                        [1.7 msec]
;d10: incremented delay (F2 in 3D) = d23/2-p14/2
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d15: TOCSY mixing time                                [12 msec]
;d16: delay for homospoil/gradient recovery
;d20: incremented delay (F1 in 3D)                     [3 usec]
;d21: 1/(2J(CaCO))                                     [3.6 msec]
;d22: 1/(2J'(CaCO)                                     [4.4 msec]
;d23: constant time delay T(N) = 1/(4J'(NCO)           [12.4 msec]
;d25: 1/(2J'(NH))                                      [5.5 msec]
;d26: 1/(4J(NH))                                       [2.3 msec]
;d28: decremented delay (F1 in 3D) = d4+p14+d20
;d29: incremented delay (F2 in 3D) = d23/2-p14/2-p26-d25-4u
;d30: decremented delay (F2 in 3D) = d23/2-p14/2
;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)
;cnst23: Caliphatic chemical shift (offset, in ppm)
;o2p: Caliphatic chemical shift (cnst23)
;l1: loop for DIPSI2 cycle: ((p6*115.112) * l1) = mixing time
;inf1: 1/SW(Hali) = 2 * DW(Hali)
;inf2: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(Hali)) =  DW(Hali)
;nd0: 2
;in10: 1/(2 * SW(N)) = DW(N)
;nd10: 4
;in20: = in0 - in28
;in28: = d28 * 2 / td1
;in29: = in10
;in30: = in10
;NS: 16 * n
;DS: >= 16
;td1: number of experiments in F1       td1 min = 2 * d28 / in0
;td2: number of experiments in F2       td2 max = 2 * d30 / in30
;FnMODE: States-TPPI in F1 (not suitable for TPPI)
;FnMODE: States-TPPI (or TPPI) in F2
;cpds1: decoupling according to sequence defined by cpdprg1
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd1: f1 channel - 90 degree pulse for decoupling sequence
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    gp 1 : gp 2 : gp 3 : gp 4
;                         50 :   40 :   60 :   30

;for z-only gradients:
;gpz1: 50%
;gpz2: 40%
;gpz3: 60%
;gpz4: 30%

;use gradient files:   
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100

;set pl9 to 120dB when presaturation is not required
;   use pl1 + 75 to 80dB to reduce radiation damping



;Processing

;SR(F1): 1/4 SWH(F1)



;$Id: hccconhgpwg3d2,v 1.7.2.2 2008/11/06 17:04:05 ber Exp $
