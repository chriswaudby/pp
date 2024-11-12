; 13C-XL-ALSOFAST-HMQC
;
;
; P. Roessler, D. Mathieu, A. D. Gossert, Apr. 2019
;
;2D H-1/X correlation via heteronuclear zero and double quantum
; coherence
;phase sensitive EA
;with gradient selection
;with selection of 13C-labeled moieties with editing step (ALSOFAST)
;with short INEPT and delayed decoupling during acquisition
;16 step phase cycle
;
;L. Mueller,
; J. Biomol. NMR 42, 129 - 137 (2008)
;P.Schanda and B. Brutscher,
; J. Am. Chem. Soc. 127, 8014 (2005)

prosol relations=<triple>

;--- Pulse program is written with syntax for AVNEO/TS4.x, the following definitions are for compatibility with AVIII/TS3.x
#ifdef TS3
#define START_NEXT_SCAN setrtp1|0
#define DWELL_RELEASE
#define DWELL_HOLD setrtp1^0
#endif


;--- AVNEO/TS4.x will use definitions for DWELL_RELEASE etc. from Avance.incl below
#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>
#include <De.incl>

;plw1 : power for 1H
;plw2 : power for 13C
;plw12 : power for 13C garp decoupling

;p1 : 90 degree hard pulse 1H
;p2 : 180 degree hard pulse 1H
;p3 : 90 degree hard pulse 13C
;p4 : 180 degree hard pulse 13C
;p16 : short gradient pulse (450 us)
;pcpd2 : 90 degree cpd-pulse 13C (garp4, 85us)

;d1 : relaxation delay (0.2 - 0.5 s)
;d2 : INEPT delay 1/4J (~1.7m)
;d21 : delay in first INEPT <= 1.7 ms (e.g. 1.1 ms)
;d22 : delay in first INEPT <= 1.7 ms (e.g. 0.7 ms)
;d23 : delay before decoupling
;d16 : delay for homospoil/gradient recovery

;in0 : 1/(SW) (Hz)

;use gradient ratio: gp 2 : gp 3 : gp4
; -40 : 40 : 20.1
;or -80 : 80 : 40.2

;for z-only gradients:
;gpz1: -31
;gpz2: -80%
;gpz3: 80%
;gpz4: 40.2% for C-13, 16.2% for N-15

;use gradient files:
;gpnam0: SMSQ10.100
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100

define delay INEPT1
define delay INEPT2

"p2=p1*2"
"p4=p3*2"

"in0=inf1/2"

"d0=3u"

"DELTA1=d0*2+p2"

"d24=d22-(d16+d0*2+p4/2+p3+p2+p16)"

"INEPT1=d21-(p16+10u+d16+p4/2)"
"INEPT2=d22*2-(p16+d24+d16+p4+p1+d0+p3+10u)"

"d23=d2*2-d22*2"
"d31=aq-d23-10u"

"acqt0=0"

baseopt_echo

1 100u ze
  10u pl1:f1
2 1m
  10u BLKGRAD
  10u do:f2
  d1
  50u
  50u UNBLKGRAD
  50u
  p16:gp0
  5m
  20u pl1:f1
  20u pl2:f2
;----------------------------------------first INEPT
  (p1 ph20):f1
  10u
  p16:gp1
  d16
  INEPT1
  (center(p2 ph21):f1 (p4 ph21):f2)
  INEPT1
  p16:gp1
  d16
  10u
  (p1 ph20):f1
;----------------------------------------gradient encoding on 13C and 13C evolution
  (p3 ph3):f2
  DELTA1
  p16:gp2*EA
  d24
  d16
  (p4 ph4):f2
  d0
  p16:gp3*EA
  d24
  d16
  (p2 ph5):f1
  d0
  (p3 ph4):f2
;----------------------------------------end of 13C evolution
  INEPT2
  p16:gp4
  d24
  d16 pl12:f2
  10u
;----------------------------------------acquisition
  ACQ_START(ph30,ph31)
  0.1u START_NEXT_SCAN
  0.1u REC_UNBLK
  0.05u DWELL_RELEASE
 
  d23 ;--- delay before decoupling starts
  10u cpd2:f2
  d31 ;--- d31=aq-d23-10u
  0.05u DWELL_HOLD
  0.1u REC_BLK
  rcyc=2
  1m do:f2 mc #0 to 2
    F1EA(calgrad(EA), caldel(d0, +in0) & calph(ph3, +180) & calph(ph31, +180))

  10u do:f1
  10u do:f2
  10u LOCKH_OFF
exit


ph3 =0 2
ph4 =0 0 2 2
ph5 =0 0 0 0 2 2 2 2 1 1 1 1 3 3 3 3
ph20=0
ph21=1
ph22=2
ph23=3
ph30=0
ph31=0 2 2 0 0 2 2 0 2 0 0 2 2 0 0 2