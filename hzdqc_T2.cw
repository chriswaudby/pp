;H(Z/D)QC T2
; set d20 = relaxation time
;for 13C, same processing as SFHZDQC
;run as pseudo-3D (td1 = 2), add and subtract to obtain Z/D components
;TODO check which is Z and D!
;
;Added option for off-resonance presat (e.g. to suppress urea signal), 21/6/15
;
;With option for 1D (first row)
;
;sfhmqcf2gpph
;avance-version (09/11/18)
;SOFAST HMQC
;2D H-1/X correlation via heteronuclear zero and double quantum
;   coherence
;phase sensitive
;with decoupling during acquisition
;
;P.Schanda and B. Brutscher, J. Am. Chem. Soc. 127, 8014 (2005)
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


"d11=30m"
"d12=20u"
"d13=4u"
"d21=1s/(cnst4*2)"
"d22=1s/(cnst4*8)"

"p2=p1*2"
"p4=p3*2"

"in0=inf2"
"d0=in0/2-p3*4/3.1415"
"l2=td1/2"

"l0=0"

baseopt_echo
"acqt0=0"

aqseq 312



  ze 
1 d11 pl12:f2
2 10m do:f2

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

  ; begin main sequence
  (p1 ph11):f1
  "DELTA1=d21-p1*0.6366"
  DELTA1

  (p3 ph12):f2
  if "l0 % 4 == 0" goto 11
  if "l0 % 4 == 1" goto 12
  if "l0 % 4 == 2" goto 13
  if "l0 % 4 == 3" goto 14

11  (p4 ph1):f2 ; A
    d22
    d20*0.5
    (center (p2 ph13):f1 (p4 ph1):f2)
    d22
    d20*0.5
    d0
    goto 20

12  d22         ; A'
    (p4 ph1):f2
    d20*0.5
    (p4 ph1):f2
    d22
    (p2 ph13):f1
    d20*0.5
    d0
    goto 20

13  d0           ; B
    d20*0.5
    d22
    (center (p2 ph13):f1 (p4 ph1):f2)
    d20*0.5
    d22
    (p4 ph1):f2
    goto 20

14  d0          ; B'
    d20*0.5
    (p4 ph1):f2
    d22
    (p2 ph13):f1
    d20*0.5
    (p4 ph1):f2
    d22
    goto 20

20  (p3 ph1):f2
    "DELTA2=d21-4u"
    DELTA2
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

    ; outer loop (d0)
    30u id0
    lo to 1 times l2

exit 
  

ph1=0 
ph2=0 
ph11=0 0 2 2
ph12=0 2
ph13=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
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
;NS: 2 * n
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
