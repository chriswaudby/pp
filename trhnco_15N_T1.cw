; 15N-T1 relaxation experiment with TROSY-HNCO read-out (pseudo-4D)
; based on Lakomek 2012
;
; relaxation times set as vclist * 40 ms
; use even numbers only for vclist
;
; F1(H) -> F3(N, vc [3rd dim]) -> F2(C=O,t1) -> F3(N,t2) -> F1(H,t4)
;
;on/off resonance Ca and C=O pulses using shaped pulse
;phase sensitive (t1 / 13C)
;phase sensitive using constant-time Echo/Antiecho (t2 / 15N)



#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

#define TEMP_COMPENSATION

"d11=30m"
"d23=12.5m"  ; 1/4J(NC')
"d26=2.6m"  ; 1/4J(HN)

; 13C evolution (initial zero dwell)
"in0=inf1/2"
"d0=3u"
"TAU=d0*2+larger(p14,p21*2)-p14"

; 15N evolution (initial half-dwell)
"in10=inf2/2"
"in30=inf2/2"
"d10=3u"
"d29=d23-p14/2-larger(p14,p21*2)/2-p25-d16"
"d30=d23-larger(p14,p21*2)/2-4u-p25-d16+3u+p14/2-in10-p21*4/PI"

; delays for T1 block
"d24=10m-p14-2u"
"d25=d24-0.5*p15"

; temperature compensation
"d17=d1-p18"

; delays for coherence transfer
"DELTA=d26"
"DELTA1=d26"
"DELTA2=d26-p22-p11-300u"
"DELTA3=d26-p23-p10-300u"
"DELTA4=260u-p24-p1*0.66"
"d27=p24+35u"

; loop counters for relaxation delay
"l1=0" ; current position
"l2=0" ; actual value of ncyc

"cnst18=-800"  ; temperature compensation
"cnst23=8.6"   ; Iburp2

"spoff2=0"
"spoff3=0"
"spoff5=bf2*(cnst22/1000000)-o2"
"spoff8=0"
"spoff15=bf1*(cnst23/1000000)-o1"

define list<loopcounter> ncyc=<$VCLIST>
1 ze
  1m
2 d11 do:f2
  "l2 = (trunc(ncyc[l1] + 0.1))"
  2m
  4u pl1:f1 pl2:f2 pl3:f3

  ; purge 15N after last fid
  (p21 ph0):f3
  10u

;---------temperature compensation and d1 recovery delay---------
#ifdef TEMP_COMPENSATION
  10u fq=cnst18(bf ppm):f3   ; -800 ppm
  10u pl8:f3
  (p18 ph0):f3    ; 15N pulse is applied far off-resonance
  10u
  10u fq=0:f3
  d17
#else
  d1
#endif

;------- kill steady state 15N ------------
  1m UNBLKGRAD
  10u pl3:f3

  (p21 ph0):f3
  10u
  p16:gp0
  d16

;------- first INEPT Hz-> 2HxNz -----------
  (p1 ph0):f1
  5u
  DELTA gron1 ; soft gradient to prevent radiation damping
  5u  groff
  (center(p1*2 ph0):f1 (p21*2 ph0):f3)
  5u
  DELTA gron1
  5u groff

;------- rephase  2HxNz to Nz------ --------
  (p1 ph1):f1  (p21 ph0):f3
  5u
  DELTA1 gron2 ; soft gradient to prevent radiation damping
  5u groff
  (center (p1*2 ph0):f1 (p21*2 ph0):f3)
  5u
  DELTA1 gron2
  5u groff
  (p21 ph11):f3 ; phase-cycle Nz, -Nz for Freeman-Hill decay
  5u
;--------------------------------------------
  (p1 ph2):f1 ; purge pulse to kill any residual HzNz
  5u
  p16:gp3  ; cleaning gradient
  d16
;------15N T1 relaxation period--------------

;if "l2==1" goto 6  ; jump to 77 for first relaxation data point, needs to be 0 in vclist
5 d24
  (p14:sp3 ph0 4u p14:sp5 ph0):f2  ; 180 C' on-res / 180 CA off-res
  d25
  (p15:sp15 ph0):f1
  d25
  (p14:sp3 ph0 4u p14:sp5 ph0):f2  ; 180 C' on-res / 180 CA off-res
  d25
  lo to 5 times c   ; delay=c*2*d25 (20ms)

;------transfer Nz to NzCz -------------------
6 3u
  (p21 ph0):f3
  d23
  (center (p14:sp3 ph0):f2 (p21*2 ph0):f3 )  ; 180 C' on-res
  d23
  (p21 ph5):f3

  p16:gp4  ; cleaning gradient
  d16

;------C' evolution (t1) ---------------------
;  (p13:sp2 ph12):f2  ; 90 C' on-res
;  d0
;  (center (p14:sp5 ph0):f2 (p21*2 ph0):f3)  ; 180 CA off-res
;  d0
;  (p13:sp8 ph0):f2  ; 90 TR C' on-res
  (p13:sp2 ph12):f2  ; 90 C' on-res
  d0
  (center (p14:sp5 ph0):f2 (p21*2 ph0):f3)  ; 180 CA off-res
  d0
  4u
  (p14:sp3 ph0):f2
  TAU
  (p14:sp5 ph0):f2
  4u
  (p13:sp8 ph0):f2

;------transfer NzCz back to Nz with constant-time 15N evolution (t2) ----------
  (p21 ph14):f3
  d10
  (p14:sp5 ph0):f2
  d29
  p25:gp5*EA
  d16
  (center (p14:sp3 ph0):f2 (2*p21 ph0):f3 )
  4u
  p25:gp5*EA*-1
  d16
  d30

;------ start TROSY read-out------------------------------------
  (p1 ph1):f1    ; Echo
  3u
  3u pl0:f1
  (p11:sp11 ph11:r):f1
  6u
  5u pl1:f1
;goto 9 ; optimization of  water supression
  DELTA2
  p22:gp2
  300u
  (center (p1*2 ph0):f1 (p21*2 ph0):f3)
  7u
  p22:gp2
  DELTA2
  300u pl0:f1
;-------------------------------------------------
  (p11:sp12 ph2):f1
  5u
  3u pl1:f1
  (p1 ph0):f1 (p21 ph8):f3  ; Echo
;goto 9 ; for optimization of  water supression
  DELTA3
  p23:gp8
  200u
  100u pl10:f1
  (center(p10 ph10:r 5u pl1 p1*2 ph0 5u pl10 p10 ph10:r):f1 (p21*2 ph0 d27):f3)
  5u
;goto 9 ; for optimization of  water supression
  p23:gp8
  DELTA3
  DELTA4
  (p21 ph0):f3
  5u
  p24:gp6 ; Echo/Anti-echo decoding gradient
9 5u
  5u pl31:f2
  20u BLKGRAD
  go=2 ph31 cpd2:f2
  d11 do:f2 mc #0 to 2
     F1PH(calph(ph12, +90), caldel(d0, +in0))
     F2EA(calgrad(EA) & calph(ph6, +180) & calph(ph7, +180) & calph(ph8, +180), caldel(d10, +in10) & caldel(d29, +in29) & caldel(d30, -in30))
     F3QF(calclc(l1, 1))

1m
1m BLKGRAD
exit

ph0=0
ph1=1
ph2=2
ph3=3
ph5=1 ; 3?
ph6=1
ph7=3
ph8=1
ph10=2
ph11=1 1 3 3
ph12=0 2
ph13=1
ph31=1 3 3 1


;-------------NOTES----------------------

;o1p = 4.7 ppm
;o2p=176 ppm (CO)
;o3p=119 ppm

;NS=8*n
;in0=inf/2
;SW=1/(2*in0)
;echo-antiecho in N15 (process as Complex in NmrDraw before splitting the spectra)

; loop counters
;l3: number of complex points (td1 / 2)
;l6: number of relaxation points

;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                               [30 msec]
;d16: delay for homospoil/gradient recovery
;d17: relaxation delay (subtracting temperature compensation pulse)

; 1H pulses

;p1:  90 deg hard 1H pulse @pl1
;pl1: 1H 90 deg
;pl0: 120 dB
;p10: 1000u (@ 700 MHz) 90 deg soft rectangular water flip-back pulse (pl10)
;p11: 1600u (@ 700 MHz) 90 deg Sinc1.1000 water flip-back pulse (sp11,sp12)
;p15: 1700u (@ 700 MHz) 180 deg IBurp2 pulse on 1H (sp15)
;sp15: 180 deg IBurp2 pulse on 1H (p15)
;sp11: 90 deg Sinc1.1000 water flip-back pulse
;sp12: 90 deg Sinc1.1000 water flip-back pulse
;spnam15: IBurp2
;spnam11: Sinc1.1000
;spnam12: Sinc1.1000
;spoffs15: 2730Hz @ 700 MHz (8.6 ppm) , should be centered in amide region but not touch the water
;cnst23: 8.6 ppm offset for Iburp2
; 13C pulses

;p13: f2 channel -  90 degree shaped pulse
;p14: f2 channel - 180 degree shaped pulse
;sp2: f2 channel - shaped pulse  90 degree  (C=O on resonance)
;sp3: f2 channel - shaped pulse 180 degree  (C=O on resonance)
;sp5: f2 channel - shaped pulse 180 degree  (Ca off resonance)
;sp8: f2 channel - shaped pulse  90 degree  (C=O on resonance)
;                  for time reversed pulse

;cnst22: Calpha chemical shift (offset, in ppm)[56 ppm]

;CPDPRG2: garp (aq C' decoupling)
;pcpd5: C' decoupling (140u or 280u @pl31)
;pl31: C' decoupling power

;15N pulses
;p21: 90 deg hard 15N pulse @pl7
;p18 : maximum duration of spin-lock; temperature compensation
;pl3 :15N 90 deg
;pl8: 15N spin-lock power

; gradients
;p16: homospoil/gradient pulse                         [1 msec]
;p22: 300u
;p23: 1000u
;p24: 60.8u Echo/Anti-echo decoding gradient
;p25: 300u Echo/Anti-echo half-encoding gradient

;for z-only gradients
;gpz0: 53%
;gpz1: 3%
;gpz2: 2%
;gpz3: 43%
;gpz4: 31%
;gpz5: -33%
;gpz6: 33%
;gpz7: -10%
;gpz8: 17%

;gpnam0 SMSQ10.100
;gpnam3 SMSQ10.100
;gpnam4 SMSQ10.100
;gpnam5 SINE.10
;gpnam6 SINE.10
;gpnam7 SMSQ10.32
;gpnam8 SMSQ10.100
