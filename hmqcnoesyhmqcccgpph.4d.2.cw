;4D 13C HMQC-NOESY-13C HMQC
;gradient selection in second HMQC
;for methyl-methyl NOES
;NUS with random quadrature detection
; (thanks to Daniel Nietlispach)
;Chris W, Dec 2017
;
; uses vclist and vplist to set evolution times and phases
; vclist = t1
;          t2
;          t3
;          ...
; vplist = 1 (r)
;          2 (i)
;          1 (E/AE)
;          ...
;
;F1(H) -> F2(H[mq],t1,d30) -> F2(C[mq],t2,d0) ---NOE--> F1(H) -> F2(C[mq],t3,d10) -> F1(H,t4)
;
;MQ evolution for 1H (taking advantage of methyl trosy)
;Uses half-dwell first-point delay by default in all indirect dims
;Option for off-res presat
;Removal of 13C equilibrium magnetisation
;Delays adjusted for zero first-order phase correction in acqusition dim


;$CLASS=HighRes
;$DIM=4D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

;------------standard delays/pulses
"p2=p1*2"
"p4=p3*2"
"d2=1s/(cnst2*2)"
"d11=30m"


;------------summary
;               t1      t2      t3      t4
;nucleus        1H      13C     13C     1H
;               noe     noe     dir     dir
;axis           F1      F2      F3      F4
;delay          d0      d10     d20
;increments     id0     id10    id20
;initial dwell  1/2     1/2     1/2
;mode           states  states  E/AE    direct
;tppi?          tppi    tppi    tppi
;aqseq 4321

;------------increments
"in0=inf1/2"  ; t1 (1Hnoe)
"in10=inf2/2"  ; t2 (13Cnoe)
"in20=inf3/2"  ; t3 (13Cdir)

;------------initial evolution times
"d0=in0/2-p3"  ; t1 (1Hnoe)
"d10=in10/2-p2"  ; t2 (13Cnoe)
"d20=3u"  ; t3 (13Cdir)

;------------other delays
"TAU=d8-p16*2-d16*2-p3"
"DELTA=p3*0.6366+p17+d17-p2-d0"
"DELTA1=d2-p16-d16"
"DELTA2=d2-0.6366*p1"
"DELTA3=d2-p17-d16-4u-de"

"acqt0=de"


;------------NUS preamble
; loop counters
"l5=1" ; t1
"l6=1" ; t2
"l7=1" ; t3
"l0=1" ; for gradient selection in t3

;set l1 = numner of NUS points to record



;------------loop for NUS points
ze
1 d11 pl12:f2

;------------calculate NUS stuff here
; 1. reset loop counters and phases
20u rp11
20u rp12
20u rp13
20u rp31
20u ru5
20u ru6
20u ru7

; 2. set l5,l6,l7 = vc
41 20u iu5
  lo to 41 times c
  20u ivc

51 20u iu6
  lo to 51 times c
  20u ivc

61 20u iu7
  lo to 61 times c
  20u ivc

"d31=20m-20u*(l5-1)-20u*(l6-1)-20u*(l7-1)"
d31 ; compensate for calculation time

; 3. states-tppi for phases p11,p12,p13
if "l5%2 == 1"
  {
  ip11*2
  ip31*2
  }
if "l6%2 == 1"
  {
  ip12*2
  ip31*2
  }
if "l7%2 == 1"
  {
  ip13*2
  ip31*2
  }

; 4. set phases for quadrature detection
if "vp > 1"  ; t1
  {
  ip11
  }
  20u ivp

if "vp > 1"  ; t2
  {
  20u ip12
  }
  20u ivp

if "vp == 1"  ; t3 (echo/anti-echo)
  {
  "l0=1"
  }
else
  {
  "l0=2"
  }
  20u ivp

; 5. calculate evolution times
20u rd0
20u rd10
20u rd20
"d0=d0+in0*(l5-2)"
"d10=d10+in10*(l6-2)"
"d20=d20+in20*(l7-2)"
20u



;------------ns loop
2 d11 do:f2
  4u BLKGRAD

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/
  4u pl9:f1
  d1 cw:f1 ph29
  4u do:f1
  4u pl1:f1 pl2:f2
  4u fq=0:f1
  4u UNBLKGRAD

;------------kill equm 13C magnetisation
  (p3 ph1):f2
  4u
  p16:gp1
  d16*2

;------------start first 13C HMQC element
  (p1 ph11):f1
  DELTA1
  p16:gp2
  d16

  (p3 ph12):f2
  ; t1 - 1H F1 indirect evolution (as MQ)
  d0
  (p4 ph2):f2
  d0
  ; t2 - 13C F2 evolution (MQ)
  d10
  (p1 ph1):f1
  (p2 ph2):f1
  (p1 ph1):f1
  d10
  (p3 ph1):f2

  p16:gp2
  d16
  DELTA1
  (p1 ph3):f1


;------------NOE mixing period
  TAU

  p16:gp3*0.71
  d16
  (p3 ph1):f2
  p16:gp3
  d16


;------------start second 13C HMQC element (gradient selected)
  (p1 ph1):f1
  DELTA2

  (p3 ph13):f2
  p17:gp4
  d17
  (p4 ph2):f2
  DELTA
  d20
  (p1 ph1):f1
  (p2 ph2):f1
  (p1 ph1):f1
  d20
  DELTA
  (p4 ph1):f2
  p17:gp4
  d17
  (p3 ph4):f2

  DELTA3 pl12:f2
  if "l0 % 2 == 1"
    {
    p17:gp5
    d16
    }
  else
    {
    p17:gp5*-1
    d16
    }
  4u BLKGRAD

;------------acquisition (ns)
  go=2 ph31 cpd2:f2

;------------loop back and calculate next NUS/RQC point
  d11 do:f2 wr #0 if #0 zd
  lo to 1 times l1

  20u BLKGRAD

exit


ph1= 0
ph2= 1
;ph3=(8) 1
ph3= 0
ph11=0
ph12=0 2
ph13=0
ph29=0
ph31=0 2


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse                       [1 msec]
;p17: second gradient pulse                          [250 usec]
;d0 : incremented delay (4D)
;d10: incremented delay (4D)
;d20: incremented delay (4D)
;d1 : relaxation delay; 1-5 * T1
;d2 : 1/(2J)CH
;d8 : mixing time
;d11: delay for disk I/O                             [30 msec]
;d16: delay for homospoil/gradient recovery [200 usec]
;d17: short delay for homospoil/gradient recovery [100 usec]
;l1: number of NUS points
;cnst2: = J(CH)
;cnst21: frequency in Hz for off-res presat
;inf1: 1/SW(H)
;inf2: 1/SW(C) = 2 * DW(C)
;inf3: 1/SW(C) = 2 * DW(C)
;in0: 1/(2 * SW(H)) = DW(H)
;in10: 1/(2 * SW(C)) = DW(C)
;in20: 1/(2 * SW(C)) = DW(C)
;nd0: 2
;nd10: 2
;nd20: 2
;NS: 2 * n
;DS: 16
;td1: number of experiments in F1
;td2: number of experiments in F2
;td3: number of experiments in F3
;FnMODE: States-TPPI (or States) in F1
;FnMODE: States-TPPI (or States) in F2
;FnMODE: States-TPPI (or States) in F3
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 31%
;gpz2: 13%
;gpz3: 50%
;gpz4: 40%
;gpz5: 20.1%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SINE.20
;gpnam5: SINE.20

                                          ;preprocessor-flags-start
;OFFRES_PRESAT: for off-resonance presaturation, set cnst21=o1(water)
;F2_plane: for zero 13C phase evolution in F3
;F3_plane: for zero 13C phase evolution in F2
;NUS: for non-uniform sampling (Topspin 3)
                                          ;preprocessor-flags-end
