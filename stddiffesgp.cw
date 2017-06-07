;Using trip1e proso1 re1ations (square pulses for excitation sculpting)
;Power level for square excitation-sculpting pulses hard-coded
;Can apply fine adjustment using cnst63
;use baseopt (Chris W Oct 2016)
;using F1 channel
;
;Using pl10 for TOCSY spin-lock (rather than pl29)
;pl10 also power-level for trim pulses
;
;stddiffesgp.3
;avance-version (10/02/02)
;pseudo 2D sequence
;   for saturation transfer difference
;with shaped pulse train for saturation on f2 channel
;alternating between on and off resonance
;   to be defined by fq2list
;with spoil sequence to destoy unwanted magnetization
;water suppression using excitation sculpting with gradients
;with spinlock to suppress protein signals
;(use parameterset STDDIFFESGP.3)
;
;M. Mayer & B. Meyer, Angew. Chem. Int. Ed. 38, 1784-1788 (1999)
;M. Mayer & B. Meyer, Angew. Chem. 111, 1902-1906 (1999).
;T.-L. Hwang & A.J. Shaka, J. Magn. Reson.,
;   Series A 112 275-279 (1995)
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


define list<frequency> stdlist=<$FQ1LIST>


"p2=p1*2"
"d12=20u"
"d11=30m"


"p29=d29"

"l5=d20/p13+0.5"
"d31=p13*l5"

"DELTA1=d1-d31"

"TAU=p1*2/3.1416+50u"
"acqt0=0"


1 ze
  3m st0
2 6m
3 6m
4 d11
  6m

5 50u UNBLKGRAD
  4u pl10:f1
  (p17 ph2)
  (p17*2 ph3)
  4u
  p30:gp1
  10m pl1:f1
  4u BLKGRAD

  DELTA1

10u  stdlist:f1

6 (p13:sp13 ph4):f1
  4u
  lo to 6 times l5
10u fq=0:f1

  p1 ph1
  4u pl10:f1
  (p29 ph5)

  50u UNBLKGRAD
  p16:gp2
  d16 pl0:f1
  (p12:sp11 ph6:r):f1
  4u
  d12 pl1:f1

  p2 ph7

  4u
  p16:gp2
  d16 
  TAU
  p16:gp3
  d16 pl0:f1
  (p12:sp11 ph8:r):f1
  4u
  d12 pl1:f1

  p2 ph9

  p16:gp3
  d16
  4u BLKGRAD

  goscnp ph31

  3m stdlist.inc
  st
  lo to 3 times nbl
  3m ipp1 ipp5 ipp6 ipp7 ipp8 ipp9 ipp31
  3m stdlist.res
  lo to 4 times ns
  d11 wr #0 
  3m rpp1 rpp5 rpp6 rpp7 rpp8 rpp9 rpp31
  3m zd
  lo to 5 times l4
exit


ph1=0 2
ph2=0
ph3=1
ph4=0
ph5=1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3
ph6=0 0 1 1
ph7=2 2 3 3
ph8=0 0 0 0 1 1 1 1
ph9=2 2 2 2 3 3 3 3
ph31=0 2 2 0 2 0 0 2


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)       [120 dB]
;pl10: f1 channel - power level for TOCSY-spinlock & trim pulse
;sp1 : f1 channel - shaped pulse 180 degree
;sp13: f2 channel - shaped pulse  for saturation          [40 - 60 dB]
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p12: f1 channel - 180 degree shaped pulse (Squa100.1000) [2 msec]
;p13: f2 channel -  shaped pulse for saturation           [50 msec]
;p16: homospoil/gradient pulse
;p17: f1 channel - trim pulse                             [2.5 msec]
;p29: f1 channel - spin-lock pulse
;p30: gradient pulse                                      [3 msec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                                  [30 msec]
;d12: delay for power switching                           [20 usec]
;d16: delay for homospoil/gradient recovery
;d20: saturation time
;d29: spinlock time                                       [10 - 50 msec]
;d31: saturation time as executed
;l4: l4 = number of averages = (total number of scans) / NS
;l5: loop for saturation: p13 * l5 = saturation time
;NS: 8 * n
;DS: 4
;td1: number of experiments
;NBL: NBL = number of irradiation frequencies

;define FQ2LIST (irradiation frequencies)
;               (list has to be stored in "/u/exp/stan/nmr/lists/f1")


;use gradient ratio:    gp 1 : gp 2 : gp 3
;                         40 :   31 :   11

;for z-only gradients:
;gpz1: 40%
;gpz2: 31%
;gpz3: 11%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100


;this pulse program produces a ser-file (PARMOD = 2D)



;The STD experiment is protected by international patents owned by:
;Alepharma Licensing, Raamfeld 67, 22397 Hamburg, Germany.
;For commercial use (direct or indirect) please contact the company for
;licensing information at:
;E-mail: info@alepharma-licensing.com,
;Fax: +49 4060847812,
;Tel: +49 1701685158 or +49 1712788867.



;$Id: stddiffesgp.3,v 1.5.2.2 2010/02/03 08:54:30 ber Exp $
