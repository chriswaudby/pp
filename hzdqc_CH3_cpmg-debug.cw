;H(Z/D)QC CPMG
; set d20 = relaxation time
; set td = 4x desired number of (real) points
;          (2x multiplet suppression, 2x ZQ/DQ selection)
;TODO check which is Z and D!
;
;Option for off-resonance presat (-DOFFRES_PRESAT)

prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"d11=30m"
"d12=20u"
"d13=4u"
"d21=1s/(cnst4*2)"
"d22=1s/(cnst4*8)"

"p2=p1*2"
"p4=p3*2"

"in0=inf2"
"d0=in0/2-p3*4/3.1415"

"l0=0"
"l1=td1"
"l2=td2/8"

define list<loopcounter> ncyc=<$VCLIST>
;"DELTA=d20"

baseopt_echo
"acqt0=0"

aqseq 312



  ze 
1 d11 pl12:f2
2 100u do:f2

  20u pl11:f1
  (2mp ph1):f1
  20u
  (3mp ph2):f1

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/

  d12 pl9:f1
  d1 cw:f1 ph29
  20u do:f1
  30u fq=0:f1

  d12 pl1:f1 pl2:f2
  50u UNBLKGRAD

  ; purge equilibrium magnetisation
  (p3 ph1):f2
  p16:gp2
  d16

  ; calculate CPMG delays
    if "ncyc == 0" goto 77
    "DELTA = d20/(4*ncyc)-larger(p2,p4)"
    ;goto 78
77  4u; "DELTA = d20*1"
78  4u

  ; reset phase cycling for XY16 CPMG
  d12 rpp21 rpp22 rpp23 rpp24

  ; begin main sequence
  (p1 ph11):f1
  "DELTA1=d21-p1*0.6366"
  DELTA1

  (p3 ph12):f2
  if "ncyc == 0" goto 20

  if "l0 % 4 == 0" goto 11
  if "l0 % 4 == 1" goto 12
  if "l0 % 4 == 2" goto 13
  if "l0 % 4 == 3" goto 14

; 2u rpp21 rpp22 rpp23 rpp24 rpp25
; first half of CPMG
;   DELTA
;   (center (p1 ph24 p2 ph21 p1 ph24):f1 (p3 ph24 p4 ph21 p3 ph24):f2 )
;   DELTA ipp21 ipp22 ipp23 ipp24 ipp25
; second half of CPMG
;   DELTA dpp21 dpp22 dpp23 dpp24 dpp25
;   (center (p1 ph23 p2 ph22 p1 ph23):f1 (p3 ph23 p4 ph22 p3 ph23):f2 )
;   DELTA

11  (p4 ph1):f2 ; A
    d22
    ; begin CPMG (A)
31  DELTA    
    (center (p1 ph1 p2 ph2 p1 ph1):f1 (p3 ph1 p4 ph2 p3 ph1):f2 )
    DELTA
    DELTA
    (center (p1 ph2 p2 ph3 p1 ph2):f1 (p3 ph2 p4 ph3 p3 ph2):f2 )
    DELTA
    lo to 31 times ncyc
; 31  DELTA
;     (center (p1 ph24 p2 ph21 p1 ph24):f1 (p3 ph24 p4 ph21 p3 ph24):f2 )
;     DELTA ipp21 ipp22 ipp23 ipp24 ipp25
;     lo to 31 times ncyc
; 41  DELTA dpp21 dpp22 dpp23 dpp24 dpp25
;     (center (p1 ph23 p2 ph22 p1 ph23):f1 (p3 ph23 p4 ph22 p3 ph23):f2 )
;     DELTA
;     lo to 41 times ncyc
; 31  DELTA
;     (center (p2 ph2):f1 (p4 ph2):f2)
;     DELTA
;     lo to 31 times 2
    ; end CPMG (A)
    (center (p2 ph13):f1 (p4 ph1):f2)
    d22
    d0
    goto 99

12  d22         ; A'
    (p4 ph1):f2
    ; begin CPMG (A')
32  DELTA    
    (center (p1 ph1 p2 ph2 p1 ph1):f1 (p3 ph1 p4 ph2 p3 ph1):f2 )
    DELTA
    DELTA
    (center (p1 ph2 p2 ph3 p1 ph2):f1 (p3 ph2 p4 ph3 p3 ph2):f2 )
    DELTA
    lo to 32 times ncyc
; 32  DELTA
;     (center (p1 ph24 p2 ph21 p1 ph24):f1 (p3 ph24 p4 ph21 p3 ph24):f2 )
;     DELTA ipp21 ipp22 ipp23 ipp24 ipp25
;     lo to 32 times ncyc
; 42  DELTA dpp21 dpp22 dpp23 dpp24 dpp25
;     (center (p1 ph23 p2 ph22 p1 ph23):f1 (p3 ph23 p4 ph22 p3 ph23):f2 )
;     DELTA
;     lo to 42 times ncyc
; 32  DELTA
;     (center (p2 ph2):f1 (p4 ph2):f2)
;     DELTA
;     lo to 32 times 2
    ; end CPMG (A')
    (p2 ph13):f1
    d22
    (p4 ph1):f2
    d0
    goto 99

13  d0           ; B
    ; begin CPMG (B)
33  DELTA    
    (center (p1 ph1 p2 ph2 p1 ph1):f1 (p3 ph1 p4 ph2 p3 ph1):f2 )
    DELTA
    DELTA
    (center (p1 ph2 p2 ph3 p1 ph2):f1 (p3 ph2 p4 ph3 p3 ph2):f2 )
    DELTA
    lo to 33 times ncyc
; 33  DELTA
;     (center (p1 ph24 p2 ph21 p1 ph24):f1 (p3 ph24 p4 ph21 p3 ph24):f2 )
;     DELTA ipp21 ipp22 ipp23 ipp24 ipp25
;     lo to 33 times ncyc
; 43  DELTA dpp21 dpp22 dpp23 dpp24 dpp25
;     (center (p1 ph23 p2 ph22 p1 ph23):f1 (p3 ph23 p4 ph22 p3 ph23):f2 )
;     DELTA
;     lo to 43 times ncyc
; 33  DELTA
;     (center (p2 ph2):f1 (p4 ph2):f2)
;     DELTA
;     lo to 33 times 2
    ; end CPMG (B)
    d22
    (center (p2 ph13):f1 (p4 ph1):f2)
    d22
    (p4 ph1):f2
    goto 99

14  d0          ; B'
    ; begin CPMG (B')
34  DELTA    
    (center (p1 ph1 p2 ph2 p1 ph1):f1 (p3 ph1 p4 ph2 p3 ph1):f2 )
    DELTA
    DELTA
    (center (p1 ph2 p2 ph3 p1 ph2):f1 (p3 ph2 p4 ph3 p3 ph2):f2 )
    DELTA
    lo to 34 times ncyc
; 34  DELTA
;     (center (p1 ph24 p2 ph21 p1 ph24):f1 (p3 ph24 p4 ph21 p3 ph24):f2 )
;     DELTA ipp21 ipp22 ipp23 ipp24 ipp25
;     lo to 34 times ncyc
; 44  DELTA dpp21 dpp22 dpp23 dpp24 dpp25
;     (center (p1 ph23 p2 ph22 p1 ph23):f1 (p3 ph23 p4 ph22 p3 ph23):f2 )
;     DELTA
;     lo to 44 times ncyc
; 34  DELTA
;     (center (p2 ph2):f1 (p4 ph2):f2)
;     DELTA
;     lo to 34 times 2
    ; end CPMG (B')
    (p4 ph1):f2
    d22
    (p2 ph13):f1
    (p4 ph1):f2
    d22
    goto 99

    ; ncyc = 0 - no relaxation delay
20  0.1u
    if "l0 % 4 == 0" goto 21
    if "l0 % 4 == 1" goto 22
    if "l0 % 4 == 2" goto 23
    if "l0 % 4 == 3" goto 24

21  (p4 ph1):f2 ; A
    d22
    (center (p2 ph13):f1 (p4 ph1):f2)
    d22
    d0
    goto 99

22  d22         ; A'
    (p4 ph1):f2
    (p2 ph13):f1
    d22
    (p4 ph1):f2
    d0
    goto 99

23  d0           ; B
    d22
    (center (p2 ph13):f1 (p4 ph1):f2)
    d22
    (p4 ph1):f2
    goto 99

24  d0          ; B'
    (p4 ph1):f2
    d22
    (p2 ph13):f1
    (p4 ph1):f2
    d22
    goto 99


99  (p3 ph1):f2
    "DELTA2=d21-4u"
    DELTA2 pl12:f2
    4u BLKGRAD

    go=2 ph31 cpd2:f2 
    20u do:f2

    d11 wr #0 if #0
    30u zd
    
    ; inner loop (pulse positions)
    30u iu0
    lo to 1 times 4

    ; middle loop (phases)
    30u ip12
    lo to 1 times 2

    ; loop over relaxation delays
    30u ncyc.inc
    lo to 1 times l1

    ; outer loop (d0)
    30u id0
    lo to 1 times l2

exit 
  
ph1=0
ph2=1
ph3=2
ph4=3
ph11=0 0 2 2
ph12=0 2
ph13=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph20=1
ph21={0 1 0 1 1 0 1 0}^2
ph22={0 3 0 3 3 0 3 0}^2
ph23=ph12+ph20
ph24=ph11-ph20
ph29=0
ph31=0 2 2 0 2 0 0 2


;pl3 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling (low power)
;p16: homospoil/gradient pulse                       [1 msec]
;p3: f2 channel - 90 degree high power pulse
;p1: f1 channel - 90 degree
;p2: f1 channel - 180 degree 
;d0 : incremented delay (2D) = in0/2-p3*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(2J)CH
;d22 : 1/(8J)CH
;cnst4: = J(CH)
;inf1: 1/SW(C) = 2 * DW(C)
;in0: 1/ SW(C) = 2 * DW(C)
;nd0: 1
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEC
;cpd3: decoupling according to sequence defined by cpdprg3: garp4.p62


;use gradient ratio:	gp 1 : gp 2
;			  11 :    7


;for z-only gradients:
;gpz1: 11%
;gpz2:  7%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100




;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1
