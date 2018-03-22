;With additional proton 90deg pulse prior to first CPD element
;Water returned to z for all increments
;Using full States-TPPI (set spectral width accordingly)
;
;hbhaconhgp3d
;avance-version (07/04/04)
;HBHACONH
;3D sequence with
;   inverse correlation for triple resonance using inept transfer steps 
;
;      F1(H,t1) -> F2(Caliph. -> Ca) -> F2(C=O) -> F3(N,t2) -> F1(H,t3)
;
;on/off resonance Ca and C=O pulses using shaped pulse
;phase sensitive (t1)
;phase sensitive using Echo/Antiecho gradient selection (t2)
;using semi-constant time in t1
;using constant time in t2
;(use parameterset HBHACONHGP3D)
;
;(S. Grzesiek & A. Bax, J. Biomol. NMR 3, 185-204 (1993))
;(D.R. Muhandiram & L.E. Kay, J. Magn. Reson. B 103, 203-216 (1994))
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

"d3=2.2m"	; 3/(10*J_CH)
"d4=1.8m"	; 1/(4*J_CH)
"d21=12.4m"
"d22=3.6m"
"d23=3.6m"
"d24=4.4m"
"d25=5.5m"
"d26=2.3m"
"d27=12.4m"


"d0=d4"
"d10=d21/2-p14/2"
"d20=3u"
"d28=d4+p14+d20"
"d29=d21/2-p14/2-p26-d25-4u"
"d30=d21/2-p14/2"

"in0=inf1/2"
"in10=inf2/4"

"FACTOR1=trunc(d28*10000000*2/(td1-2))"
"in28=FACTOR1/10000000"

"in20=in0-in28"
"in29=in10"
"in30=in10"


"DELTA1=d22-d3"
"DELTA2=d23-p14"
"DELTA3=d27-d24+4u"
"DELTA4=d25-p16-d16-4u"
"DELTA5=p16+d16+8u"


"spoff2=0"
"spoff3=0"
"spoff5=bf2*((cnst21-cnst23)/1000000)"
"spoff7=bf2*((cnst22-cnst21)/1000000)"
"spoff8=0"


aqseq 321


1 ze
  d11 pl16:f3
2 d11 do:f3
3 d11 fq=cnst23(bf ppm):f2
  d1
  50u UNBLKGRAD
  d12 pl1:f1 pl0:f2 pl3:f3

  (p1 ph3):f1
  d0
  (p14:sp3 ph1):f2
  d20
  (p2 ph1):f1 
  d28
  (p1 ph2):f1

  (p13:sp2 ph1):f2
  d3 pl19:f1
  (p26 ph10):f1
  DELTA1 cpds1:f1 ph2
  (p14:sp3 ph1):f2
  d22
  (p13:sp8 ph1):f2

  4u
  (p14:sp5 ph1):f2
  DELTA2
  (p14:sp3 ph1):f2
  4u
  (p14:sp5 ph1):f2
  DELTA2
  (p13:sp2 ph1):f2

  4u do:f1
  (p26 ph9):f1
  p16:gp1
  d16 fq=cnst21(bf ppm):f2
  (p26 ph2):f1
  20u cpds1:f1 ph1

  (p13:sp2 ph4):f2
  d24
  (p14:sp7 ph1):f2
  DELTA3
  (center (p14:sp3 ph1):f2 (p22 ph1):f3 )
  d27
  (p14:sp7 ph1):f2
  4u
  (p13:sp8 ph1):f2

  (p21 ph1):f3
  d30
  (p14:sp7 ph1):f2
  d30
  (center (p14:sp3 ph1):f2 (p22 ph8):f3 )
  d10
  (p14:sp7 ph1):f2
  d29
  4u do:f1
  (p26 ph7):f1
  4u
  p16:gp2*EA
  d16 
  DELTA4 pl1:f1

  (center (p1 ph1):f1 (p21 ph5):f3 )
  d26
  (center (p2 ph1):f1 (p22 ph1):f3 )
  d26
  (center (p1 ph2):f1 (p21 ph6):f3 )
  d26
  (center (p2 ph1):f1 (p22 ph1):f3 )
  d26
  (p1 ph1):f1

  DELTA5
  (p2 ph1):f1
  4u
  p16:gp3
  d16 pl16:f3
  4u BLKGRAD
  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
     F1PH(rd10 & rd29 & rd30 & ip3 & dp10, id0 & id20 & dd28 & rp10 & ip9*2)
     F2EA(igrad EA & ip5*2, id10 & id29 & dd30)
exit


ph1=0
ph2=1
ph3=0 2
ph4=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph5=0 0 2 2
ph6=3 3 1 1
ph7=3
ph8=0 0 0 0 2 2 2 2
ph9=0 2
ph10=1
ph31=0 2 2 0 0 2 2 0 2 0 0 2 2 0 0 2


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl16: f3 channel - power level for CPD/BB decoupling
;pl19: f1 channel - power level for CPD/BB decoupling
;sp2: f2 channel - shaped pulse  90 degree  (on resonance)
;sp3: f2 channel - shaped pulse 180 degree  (on resonance)
;sp5: f2 channel - shaped pulse 180 degree  (C=O off resonance)
;sp7: f2 channel - shaped pulse 180 degree  (Ca off resonance)
;sp8: f2 channel - shaped pulse  90 degree  (on resonance)
;                  for time reversed pulse
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p13: f2 channel -  90 degree shaped pulse
;p14: f2 channel - 180 degree shaped pulse
;p16: homospoil/gradient pulse                         [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p26: f1 channel -  90 degree pulse at pl19
;d0 : incremented delay (F1 in 3D) = d4
;d1 : relaxation delay; 1-5 * T1
;d3 : 3/(10J(CH))                                      [2.2 msec]
;d4 : 1/(4J(CH))                                       [1.8 msec]
;d10: incremented delay (F2 in 3D) =  d21/2-p14/2
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d16: delay for homospoil/gradient recovery
;d20: decremented delay (F1 in 3D)                     [3 usec]
;d21: 1/(4J(NCO)), T(N)                                [12.4 msec]
;d22: 1/(8J(CaCb))                                     [3.6 msec]
;d23: 1/(4J'(CaCO))                                    [3.2 msec]
;d24: 1/(4J(CaCO))                                     [4.4 msec]
;d25: 1/(2J(NH))                                       [5.5 msec]
;d26: 1/(4J(NH))                                       [2.3 msec]
;d27: 1/(4J'(NCO))                                     [12.4 msec]
;d28: incremented delay (F1 in 3D) = d4+p14+d20
;d29: incremented delay (F2 in 3D) = d21/2-p14/2-p26-d25-4u
;d30: decremented delay (F2 in 3D) = d21/2-p14/2
;cnst21: CO chemical shift (offset, in ppm)
;cnst22: Calpha chemical shift (offset, in ppm)
;cnst23: Caliphatic chemical shift (offset, in ppm)
;o2p: Caliphatic chemical shift (cnst23)
;inf1: 1/SW(Hali) = 2 * DW(Hali)
;inf2: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(Hali)) =  DW(Hali)
;nd0: 2
;in10: 1/(4 * SW(N)) = (1/2) DW(N)
;nd10: 4
;in20: = in0 - in28
;in28: = d28 *2 / td1
;in29: = in10
;in30: = in10
;NS: 4 * n
;DS: >= 16
;td1: number of experiments in F1       td1 min = 2 * d28 / in0
;td2: number of experiments in F2       td2 max = 2 * d30 / in30
;FnMODE: States-TPPI in F1 (not suitable for States or TPPI)
;FnMODE: echo-antiecho in F2
;cpds1: decoupling according to sequence defined by cpdprg1
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd1: f1 channel - 90 degree pulse for decoupling sequence
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;use gradient ratio:    	gp 1 : gp 2 : gp 3
;				  30 :   80 :  8.1

;for z-only gradients
;gpz1: 30%
;gpz2: 80%
;gpz3: 8.1%

;use gradient files:
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100





;$Id: hbhaconhgp3d,v 1.14 2007/04/11 13:34:29 ber Exp $
