;15N TROSY CPMG
; Vallurupalli 2007 PNAS 104 18473-77, Fig 2
; (NB phases in paper are varian!)
;
;using f3 - channel
;(use parameterset )
;
;T. Schulte-Herbrueggen & O.W. Sorensen, J. Magn. Reson. 144,
;   123 - 128 (2000)
;(M. Czisch & R. Boelens, J. Magn. Reson. 134, 158-160 (1998) )
;(K. Pervushin, G. Wider & K. Wuethrich, J. Biomol. NMR 12,
;   345-348 (1998) )
;(A. Meissner, T. Schulte-Herbrueggen, J. Briand & O.W. Sorensen, Mol. Phys. 96,
;   1137-1142 (1998) )
;(J. Weigelt, J. Am. Chem. Soc. 120, 10778-10779 (1998) )
;(M. Rance, J.P. Loria & A.G. Palmer III, J. Magn. Reson. 136, 91-101 (1999) )
;(G. Zhu, X.M. Kong & K.H. Sze, J. Biomol. NMR 13, 77-81 (1999) )
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define delay TIME_T2
"TIME_T2 = 40m"
define list<loopcounter> ncyc = {0 40 1 36 2 32 3 28 4 24 5 20 6 18 7 16 8 14 9 12 10}

"p2=p1*2"
"p22=p21*2"
"p24=p23*2"
"d21=p21"
"d11=30m"
"d12=20u"
"d13=4u"

"in0=inf2/2"
"d0=20u"

"DELTA1=2.25m-p16-d16-larger(p1,p21)"  ; tau_a
"DELTA2=2.68m-p16-d16-larger(p2,p24)" ; tau_b
"DELTA3=2.25m-0.1u-p10-p19-d16-larger(p1,p21)"  ; tau_a
"DELTA4=2.25m-0.1u-p10-p19-d16-larger(p1,p21)-4u-p21" ; tau_a
"TAU=1s" ; dummy for CPMG inter-pulse delay

"acqt0=0"
aqseq 312

"l1=0" ; loop counter for quad. detection

1 ze 
  st0 ; interleaving
2 d11
3 d11
4 d11 BLKGRAD
  d1 pl1:f1 pl3:f3

  ; calculate CPMG pulse spacing
  if "ncyc==0" goto 80
  4u
  "TAU = TIME_T2/(4*ncyc) - p23"
  4u
80 4u
  50u UNBLKGRAD

  ; purge Nz
  (p21 ph1):f3
  p16:gp0
  d16

  ; first INEPT
  (p10:sp10 ph3):f1
  4u pl1:f1
  (p1 ph21)
; goto 999 ; optimise p10/sp10
  p16:gp1
  d16
  DELTA1
  (center (p2 ph2) (p22 ph1):f3 )
  DELTA1
  p16:gp1
  d16
  (p1 ph1) 

  p16:gp2
  d16

  d8 pl23:f3 ; equilibration delay

  ; first CPMG element
  (p23 ph11):f3
  if "ncyc==0" goto 71
10 TAU
  (p24 ph3):f3
  TAU
  lo to 10 times ncyc
71 (p23 ph12):f3

  ; P element
  p16:gp3
  d16
  DELTA2
  (center (p1 ph0 p2 ph3 p1 ph0) (p23 ph21 p24 ph22 p23 ph21):f3 )
  p16:gp3
  d16
  (ralign (p1 ph0 p2 ph3 p1 ph0) (DELTA2 p23 ph13):f3)

  ; second CPMG element
  if "ncyc==0" goto 72
20 TAU
  (p24 ph3):f3
  TAU
  lo to 20 times ncyc
72 (p23 ph3):f3

  d8 pl3:f3 ; equilibration delay
#ifdef ANTITROSY
  (p1 ph0 p2 ph3 p1 ph0)
#endif
  p16:gp4
  d16

  ; t1 evolution
if "l1 % 2 == 1" goto 73
  (p3 ph14):f3
  goto 74
73 (p3 ph24):f3
74 1u
#   ifdef LABEL_CN
  d0 gron5*-1
  2u groff
  (p8:sp13 ph1):f2
  d0 gron5
  2u groff
#   else
  d0 gron5*-1
  2u groff
  d0 gron5
  2u groff
#   endif /*LABEL_CN*/
  1u

  ; back-transfer (1)
  (p1 ph15)
  0.1u
  (p10:sp12 ph23):f1
  p19:gp6
  d16
  DELTA3
  (center (p2 ph0) (p22 ph0):f3 )
  DELTA3
  p19:gp6
  d16
  (p10:sp11 ph23):f1
  0.1u pl1:f1
  (p1 ph1) (p21 ph17):f3
  ; back-transfer (2)
  p19:gp7
  d16
  4u
  d21 ; =p21
  DELTA4
  (p10:sp10 ph2):f1	; flipdown(-y), z -> -x
  0.1u pl1:f1
  (center (p2 ph2) (p22 ph8):f3 )
  (p10:sp11 ph2):f1	; flipback(-y), x -> z
  0.1u
  DELTA4
  p19:gp7
  d16 
  (p21 ph0):f3
999  4u BLKGRAD

  ; acquistion (non-interleaved)
  go=2 ph31
  d11 mc #0 to 4
     F1QF(ncyc.inc)
     F2PH(iu1 & ip5*2 & ip6*2 & ip7*2 & ip31*2, ip4*2 & ip24*2 & id0)

  ; acquisition (interleaved)
  goscnp ph31

  d11 st ncyc.inc ; cpmg innermost loop
  lo to 2 times nbl

  20u ncyc.res
  d11 ipp11 ipp12 ipp13 ipp14 ipp24 ipp31
  lo to 3 times ns ; phase cycling loop

  d11 mc #0 to 4
     F1QF()
     F2PH(iu1 & ip5*2 & ip6*2 & ip7*2 & ip31*2, ip4*2 & ip24*2 & id0)
exit 
  

ph1=0
ph2=1
ph3=2
ph4=3
ph11=0 2
ph12=3 3 1 1
ph13=0 0 2 2
ph14=3 3 1 1 2 2 0 0
ph24=3 3 1 1 0 0 2 2
ph15=3
ph17=1
#ifdef ANTITROSY
ph20=0
ph21=3
ph22=0
ph23=3
#else
ph20=2
ph21=0
ph22=1
ph23=1
#endif /* ANTITROSY */
ph31=1 3 3 1 0 2 2 0


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl23: f3 channel - power level for CPMG
;sp10: f1 channel - shaped pulse  90 degree  (flipdown)
;sp11: f1 channel - shaped pulse  90 degree  (flipback)
;sp12: f1 channel - shaped pulse  90 degree
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p10: f1 channel -  90 degree shaped pulse             [1 msec]
;p16: homospoil/gradient pulse                         [1 msec]
;p19: homospoil/gradient pulse                         [300 usec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p23: f3 channel -  90 degree pulse at CPMG power
;p24: f3 channel - 180 degree CPMG pulse
;d0 : incremented delay (2D)                           [6 usec]
;d1 : relaxation delay; 1-5 * T1
;d8 : equilibration delay [5 ms]
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(N)) = DW(N)
;nd0: 2
;NS: 8 * n
;DS: 16
;td1: number of experiments
;FnMODE: echo-antiecho


;for z-only gradients:
;gpz0: -30%
;gpz1: 10%
;gpz2: 24%
;gpz3: 4.8%
;gpz4: 10%
;gpz5: 0.5
;gpz6: 12%
;gpz7: 50%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.32
;gpnam7: SMSQ10.32


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


