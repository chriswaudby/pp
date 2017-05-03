;CA CEST, via (HA)CA(CO)NNH experiment
;
;based on (HA)CA(CO)NH, Jun 2012
;Proton-coupled in Ca-dimension
;Rephasing of Ca-Ha coupling prior to Ca evolution to give in-phase doublets
;Purging equm 13C magnetisation
;With Ca-selective refocusing pulse in rephasing element
;
;cbcaconhgp3d
;avance-version (07/06/20)
;CBCACONH
;3D sequence with
;   inverse correlation for triple resonance using inept transfer steps 
;
;      F1(H) -> F2(Caliph.,t1 -> Ca) -> F2(C=O) -> F3(N,t2) -> F1(H,t3)
;
;on/off resonance Ca and C=O pulses using shaped pulse
;phase sensitive (t1)
;phase sensitive using Echo/Antiecho-TPPI gradient selection (t2)
;using constant time in t1
;using constant time in t2
;(use parameterset CBCACONHGP3D)
;
;S. Grzesiek & A. Bax, J. Biomol. NMR 3, 185-204 (1993)
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



;list of CEST saturation frequencies
;first line in file should be zero, indicating the reference plane
;give values in Hz relative to sfo2
define list<frequency> C13sat = <$FQ1LIST>



"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"

"d4=1.7m"			;tau a [1/(4*J_CaHa)]
"d21=12.4m"			;T N
"d23=3.6m"			;T C
"d24=4.4m"			;tau d
"d25=5.5m"			;tau f
"d26=2.3m"			;tau g
"d27=12.4m"			;tau e


"d10=d21/2-p14/2"
"d29=d21/2-p14/2-p26-d25-4u"
"d30=d21/2-p14/2"

"in10=inf2/4"
"in29=in10"
"in30=in10"

"DELTA=d26-p16-d16-4u-larger(p1,p21)"
"DELTA1=d4-p16-d16-4u-p14/2"
"DELTA2=d23-p24/2+p14*2"
"DELTA3=d27-d24+4u"
"DELTA4=d25-p16-d16-4u"
"DELTA5=p16+d16+8u+de"
"DELTA6=d4-p24/2-p16-d16-4u"

"spoff2=0"
"spoff3=0"
"spoff5=bf2*((cnst21-cnst22)/1000000)"
"spoff7=bf2*((cnst22-cnst21)/1000000)"
"spoff8=0"
"spoff9=0"


; loop counter for 'reference'
"l2=1"
aqseq 312


1 ze
  d11 pl16:f3
2 d11 do:f3
3 d11 fq=cnst22(bf ppm):f2   ; set transmitter to CA

; temperature compensation
if "l2==1"
{
  4u LOCKH_ON
  ;4u fq=cnst19(bf ppm):f1  ; 1H on HA
  4u pl8:f1
  d18 cpds8:f1
  4u do:f1
  ;4u fq=0:f1
  4u LOCKH_OFF
}

  ; purge water before recycle delay
  4u pl1:f1
  4u UNBLKGRAD
  p16:gp5*0.8
  d16
  (p1 ph1):f1
  p16:gp5
  d16
  4u BLKGRAD

  ; recycle delay
  d1
  50u UNBLKGRAD
  d12 pl1:f1 pl2:f2 pl3:f3

  ; purge eq'm 13C
  (p3 ph1):f2
  p16:gp5
  d16
  (p3 ph2):f2
  p16:gp5*0.7
  d16*2 pl0:f2

  ; begin main sequence
  (p10:sp1 ph2):f1  ; WFD
  4u pl1:f1

  (p1 ph1):f1  ; HAy --> HAxCAz
  4u
  p16:gp6
  d16
  DELTA1
  (center (p2 ph1):f1 (p14:sp3 ph1):f2 )
  DELTA1
  4u
  p16:gp6
  d16
  (p1 ph2):f1
  4u
  p16:gp7   ; HAzCAz
  d16

  (p13:sp2 ph11):f2   ; HAzCAy --> CAx
  4u
  p16:gp8
  d16
  DELTA6
  (center (p2 ph1):f1 (p24:sp9 ph1):f2 )	; Ca-selective [Reburp]
  DELTA6
  4u
  p16:gp8
  d16
  (p13:sp8 ph2):f2

  4u
  p16:gp4*0.7   ; CAz cleaning gradient
  d16

  ; CEST
  if "l2==1" goto 77

  4u LOCKH_ON
  4u C13sat:f2
;  4u fq=cnst19(bf ppm):f1  ; 1H on HA
  4u pl8:f1 pl18:f2
  d18 cpds8:f1 cw:f2 ph1
  4u do:f1 do:f2
;  4u fq=0:f1
  4u fq=cnst22(bf ppm):f2
  4u LOCKH_OFF

77 4u pl1:f1 pl0:f2
  p16:gp4   ; CAz cleaning gradient
  d16

  (p13:sp2 ph1):f2   ; CAy --> CAxCOz
  DELTA2
  (p14:sp5 ph1):f2  ; CO
  4u
  (p24:sp9 ph1):f2  ; Ca-selective [Reburp]
  DELTA2
  (p14:sp5 ph1):f2  ; CO
  4u
  (p13:sp8 ph2):f2
  d12 pl19:f1

  p16:gp1   ; CAzCOz
  d16 fq=cnst21(bf ppm):f2   ; move transmitter to CO
  (p26 ph2):f1
  20u cpds1:f1 ph1

  (p13:sp2 ph1):f2
  d24
  (p14:sp7 ph1):f2
  DELTA3
  (center (p14:sp3 ph1):f2 (p22 ph1):f3 )
  d27
  (p14:sp7 ph1):f2
  4u
  (p13:sp8 ph1):f2

  (p21 ph12):f3
  d30
  (p14:sp7 ph1):f2
  d30
  (center (p14:sp3 ph1):f2 (p22 ph1):f3 )
  d10
  (p14:sp7 ph1):f2
  d29
  4u do:f1
  (p26 ph7):f1
  4u
  p16:gp2*EA
  d16 
  DELTA4 pl1:f1

  (center (p1 ph1):f1 (p21 ph13):f3 )
  4u
  p16:gp9
  d16
  DELTA
  (center (p2 ph1):f1 (p22 ph1):f3 )
  DELTA
  4u
  p16:gp9
  d16
  (center (p1 ph2):f1 (p21 ph2):f3 )
  4u
  p16:gp10
  d16
  DELTA
  (center (p2 ph1):f1 (p22 ph1):f3 )
  DELTA
  4u
  p16:gp10
  d16

  (p1 ph1):f1
  DELTA5
  (p2 ph1):f1
  4u
  p16:gp3
  d16 pl16:f3
  4u BLKGRAD

  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
     F1QF(C13sat.inc & iu2)
     F2EA(ru2 & igrad EA & ip13*2, id10 & id29 & dd30 & ip12*2 & ip31*2)
;     F1PH(rd10 & rd29 & rd30 & rp11 & rp31 & ip3, id0 & dd20) 
exit


ph1=0
ph2=1
ph7=3
ph11=0 2
ph12=0 0 2 2
ph13=0
ph31=0 2 2 0 


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
;sp9: f2 channel - Ca-selective shaped pulse 180 degree (Reburp.1000)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p13: f2 channel -  90 degree shaped pulse
;p14: f2 channel - 180 degree shaped pulse
;p16: homospoil/gradient pulse                         [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p24: f2 channel - Ca-selective 180 degree shaped pulse
;p26: f1 channel -  90 degree pulse at pl19
;d1 : relaxation delay; 1-5 * T1
;d3 : tau b : 1.1m-p14
;d4 : 1/(4J(CH)) - tau a                               [1.7 msec]
;d10: incremented delay (F2 in 3D) =  d21/2-p14/2
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d16: delay for homospoil/gradient recovery
;d21: T(N)                                             [12.4 msec]
;d22: T(C)                                             [3.6 msec]
;d23: tau c                                            [3.6 msec]
;d24: tau d                                            [4.4 msec]
;d25: tau f                                            [5.5 msec]
;d26: 1/(4J(NH)) - tau g                               [2.3 msec]
;d27: tau e                                            [12.4 msec]
;d29: incremented delay (F2 in 3D) = d21/2-p14/2-p26-d25-4u
;d30: decremented delay (F2 in 3D) = d21/2-p14/2
;cnst21: CO chemical shift (offset, in ppm)  [174 ppm]
;cnst22: Calpha chemical shift (offset, in ppm) [55 ppm]
;o2p: Calpha chemical shift (cnst22) [55 ppm]
;inf2: 1/SW(N) = 2 * DW(N)
;in10: 1/(4 * SW(N)) = (1/2) DW(N)
;nd10: 4
;in29: = in10
;in30: = in10
;NS: 4 * n
;DS: >= 16
;td2: number of experiments in F2       td2 max = 2 * d30 / in30
;FnMODE: States-TPPI (or TPPI) in F1
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
;gpz4: 17.13%
;gpz5: 47%

;use gradient files:
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam5: SINE.100



;$Id: cbcaconhgp3d,v 1.13.2.1 2007/07/04 13:41:19 ber Exp $
