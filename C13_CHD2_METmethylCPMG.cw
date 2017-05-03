;C13_CHD2_METmethylCPMG.gcm
;Version for BOTH topSpin v2 & v3
;Change calculation of l6 below to:
;	"l6=d5/d8+0.01"      ;for TopSpin v2
;	"l6=d5/d8-0.49"      ;for TopSpin v3
;
;Pulse sequence to measure 13C-methyl R2-exchange contribution to relaxation
;for isolated 13CHD2 of Methionine.
;Including optional 2H refocusing by hard 180 during long CPMG delays.
;Optional 2H decoupling during t1.
;See instructions at the end of pulse sequence
;
;       1) use vclist for different cpmg fields
;       2) interleaved acquisition ns - cpmg - hypercomplex - t1
;       3) cpmgloop = 0,2,4,6,8,...,40
;
;Based on uw_C13_CPMG_Rex_NH_fm_500_my_v5.c, which in turn is
;based on N15_CPMG_Rex_NH_fm_500.c
;
;Coded by Göran Carlström 2014-11-18
;Edited 2014-11-28 to BLKGRAD during CPMG-periods
;Copyright 2014 Göran Carlström, Lund University, Sweden
;
;Disclaimer
;The  pulse program has been tested and used on a Bruker AVANCE III HD 500 MHz
;spectrometer at Lund University, Sweden, and on a Bruker AVANCE 600 MHz
;spectrometer enquipped with a 1H,13C,15N TCI cryoprobe at AstraZeneca R&D
;Mölndal, Sweden. The sequence is made available without any expressed
;warranty. We are not liable for any potential damage that might be caused in
;connection with the pulse program. We are not responsible for any damage on
;your hardware or software caused by the use of this pulse program.
;
;   channel 1 - 1H
;   channel 2 - 13C
;   channel 4 - 2H

;$CLASS=HighRes
;$OWNER=Goran
;$DIM=pseudo_3D


prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<loopcounter> cpmglist = <$VCLIST>

aqseq 312

"in10=inf2/2"

;***** PULSES *****
"p2=2*p1"
"p4=2*p3"
"p6=2*p5"
"p17=0.5*p16"
"p18=0.9*p16"
"p19=0.3*p16"
"p26=2.0*p25"

;***** LOOPS *****
"l1=(td1)"
"l2=(td2/2)"
"l3=0"
"l6=0"
"l7=0"
"l9=1"

;***** DELAYS *****
"d4=1.0/(4.0*cnst3)"            ;taub
"d7=1.0/(4.0*cnst3)"            ;taua
"d5=d3/4.0"			;tau_CPMG
"d17=d16/2.0"
"d10=3u"
"d11=30m"
"d12=20u"
"d13=3u"
"d21=0"
"d22=0"
"d23=d5-p5"
"d30=p1*2.385+2.5*d19-p3"
"DELTA1=d7-p17-d16-d13"
"DELTA3=d4-p17-d16-d13"
;"DELTA4=d4-p18-d16-d13"
"DELTA4=d4-p18-d16-d13-d30-p3-0.6366*p1" ; CW
;"DELTA5=DELTA4-d12"
"DELTA5=d4-p18-d16-d13-d12-d30-p3"
;"DELTA6=DELTA3-p5"
"DELTA6=d4-p17-d16-6u-larger(p2,p3)-0.6366*p1" ; CW
"DELTA7=d4+p2"

"cnst31=4.0*atan(1.0)"          ;PI

"acqt0=0"

1  ze
   d11 pl12:f2
;#if defined(H2REF) || defined(H2DEC)
   d11 LOCKDEC_ON
   50u LOCKH_ON
   d11 H2_PULSE
;#endif	/* H2REF, H2DEC */

2  d11 do:f2
   d11
3  d11
4  d11
5  d12 pl1:f1		;set transmitter power
   d12 pl2:f2		;set decoupler power
;#if defined(H2REF) || defined(H2DEC)
   d12 pl4:f4		;set power for 4th channel
   d13 BLKGRAD
   d11 H2_LOCK
   6m LOCKH_OFF
;#endif	/* H2REF, H2DEC */
   d12 fq=0:f1  	;set transmitter offset
   "l3=cpmglist"
  if "l3==0" goto 6
   "d5=d3/(4.0*l3)"
   "d23=d5-p5"
;#ifdef H2REF
  if "d5<d8" goto 9
   "l6=d5/d8-0.49"	;for TopSpin v3
;   "l6=d5/d8+0.01"	;for TopSpin v2
   "l7=l6-1"
  if "d5/(2.0*l6)<p25" goto 9
   "d21=d5/(2.0*l6)-p25"
  if "d21<p5" goto 9
   "d22=d21-p5"
   "l9=0"
   goto 6
9  d13
  "l9=1"
   d13
;#endif	/* H2REF */

6  d13
;   d12 pl9:f1		;set transmitter power for presat
;   d12 fq=cnst1:f1	;set transmitter offset for presat
;  (p9 ph10):f1		;presat
;   d13
;   d12 pl1:f1		;set transmitter power
;   d12 fq=0:f1		;set transmitter offset
   d1
;#if defined(H2REF) || defined(H2DEC)
   50u LOCKH_ON
   d12 H2_PULSE
;#endif	/* H2REF, H2DEC */
   50u UNBLKGRAD
  (p3 ph10):f2		;pulse on 2nd channel
   d13
   p16:gp1
   d16
   d12 pl0:f1		;set transmitter power
  (p14:sp1 ph13):f1	;shaped pulse on transmitter
   4u
   d12 pl1:f1		;set transmitter power

;***** INEPT *****
  (p1 ph10):f1		;pulse on transmitter
   d13
   p17:gp2
   d16
   DELTA1		;delay < 1/4J(XH)
  ( center (p2 ph10):f1 (p4 ph10):f2 )
   DELTA1
   d13
   p17:gp2
   d16
  (p1 ph11):f1		;pulse on transmitter
   d13
   p16:gp5
   d17 ;pl5:f2  ;CW
   d17 BLKGRAD
;#if defined(H2REF)
   d6 pl6:f4		;delay(taueq);
;#else
;   d6
;#endif	/* H2REF */
  (p3 ph2):f2		;pulse on 2nd channel (CW)
  3u pl5:f2  ; CW
  
;***** 1st CPMG period  *****
  if "l3 == 0" goto 7
;#if defined(H2REF)
  if "l9 > 0" goto 14   ; no 2H refocusing
  if "l7 < 1" goto 13   ; short form 2H refocusing

;***** Use pulsed 2H refocusing, long form *****
11 d21
  (p26 ph10):f4		;pulse on 4th 2H channel
   d21
  lo to 11 times l7
   d21
  (p26 ph10):f4		;pulse on 4th 2H channel
   d22
  (p6 ph11):f2		;pulse on 2nd channel
   d22
  (p26 ph10):f4		;pulse on 4th 2H channel
   d21
12 d21
  (p26 ph10):f4		;pulse on 4th 2H channel
   d21
  lo to 12 times l7
  lo to 11 times l3
   goto 7

;***** Use pulsed 2H refocusing, short form *****
13 d21
  (p26 ph10):f4		;pulse on 4th 2H channel
   d22 
  (p6 ph11):f2		;pulse on 2nd channel
   d22
  (p26 ph10):f4		;pulse on 4th 2H channel
   d21
  lo to 13 times l3
   goto 7
;#endif	/* H2REF */

;***** No pulsed 2H refocusing *****
14 d23
  (p6 ph11):f2		;pulse on 2nd channel
   d23
  lo to 14 times l3

;***** change 2CyHz to Cx *****
7  d13 UNBLKGRAD
   p17:gp6
   d16 pl2:f2  ; CW
   DELTA6

;***** composite 1H 90y-180x-90y on top of 13C 180x *****
  ( center (p1 ph11 p2 ph10 p1 ph11):f1 (p4 ph10):f2 )

   DELTA6
   d13 pl5:f2
   p17:gp6
   d17  ; d17=0.5*d16
   d17 BLKGRAD

;***** start of the CPMG train for second period time_T2/2 on Cx *****
  if "l3 == 0" goto 8
;#if defined(H2REF)
  if "l9 > 0" goto 24
  if "l7 < 1" goto 23 ; CW error! originally goto 22

;***** Use pulsed 2H refocusing, long form *****
21 d21
  (p26 ph10):f4		;pulse on 4th 2H channel
   d21
  lo to 21 times l7
   d21
  (p26 ph10):f4		;pulse on 4th 2H channel
   d22
  (p6 ph10):f2		;pulse on 2nd channel
   d22
  (p26 ph10):f4		;pulse on 4th 2H channel
   d21
22 d21
  (p26 ph10):f4		;pulse on 4th 2H channel
   d21
  lo to 22 times l7
  lo to 21 times l3
   goto 8

;***** Use pulsed 2H refocusing, short form *****
23 d21
  (p26 ph10):f4		;pulse on 4th 2H channel
   d22
  (p6 ph10):f2		;pulse on 2nd channel
   d22
  (p26 ph10):f4		;pulse on 4th 2H channel
   d21
  lo to 23 times l3
   goto 8
;#endif	/* H2REF */

;***** No pulsed 2H refocusing *****
24 d23
  (p6 ph10):f2		;pulse on 2nd channel
   d23
  lo to 24 times l3

8  d13  pl2:f2 ; CW
  (p3 ph11):f2		;pulse on 2nd channel CW
   d13 UNBLKGRAD
   d6       		;delay(taueq);
   p19:gp7
;#if defined(H2DEC)
   d17 pl4:f4
   d17 BLKGRAD
;#else
;   d16
;#endif /* H2DEC */

;#if defined(H2DEC)
;***** 2D decoupling on *****
  (p23 ph11):f4		;pulse on 4th 2H channel
   4u cpd4:f4
;#endif /* H2DEC */
  (p3 ph3):f2		;pulse on 2nd channel
   d4			;delay(tau1 + taub);

;***** T1 evolution *****
   d10
  (p2 ph10):f1		;pulse on transmitter
   d10			;delay(tau1);

  (p4 ph10):f2		;pulse on 2nd channel
   DELTA7		;delay(taub + 2*pw);
  (p3 ph11):f2		;pulse on 2nd channel
;#if defined(H2DEC)
   4u do:f4
  (p23 ph13):f4		;pulse on 4th 2H channel
;***** 2D decoupling off *****
;#endif /* H2DEC */

   d13 UNBLKGRAD
   p16:gp4
   d16
  (p1 ph10):f1		;pulse on transmitter
   d13
   p18:gp3
   d16
   DELTA4
  (p1*0.231 ph11 d19 p1*0.692 ph11 d19 p1*1.462 ph11 d19):f1 (d30 p4 ph10):f2
  (p1*1.462 ph13 d19 p1*0.692 ph13 d19 p1*0.231 ph13):f1
   d13
   d12 pl12:f2		;set decoupler power
   DELTA5 
   p18:gp3
   d17
   d17 BLKGRAD
  go=2 ph31 cpd2:f2
   d11 do:f2 wr #0 if #0 zd
   d11 cpmglist.inc
  lo to 3 times l1
   d11 ip3 cpmglist.res
  lo to 4 times 2
   d11 id10
  lo to 5 times l2

;#if defined(H2REF) || defined(H2DEC)
   d11 H2_LOCK
   d11 LOCKH_OFF
   d11 LOCKDEC_OFF
;#endif	/* H2REF, H2DEC */
exit

ph2 = 0 2
ph3 = 0 0 2 2
ph4 = 0
ph10 = 0
ph11 = 1
ph12 = 2
ph13 = 3
ph31 = 0 2 2 0


;pl0 : 120dB
;pl1 : f1 channel - power level for 1H pulse (default)
;pl2 : f2 channel - power level for 13C pulse (default)
;pl4 : f4 channel - power level for 2H pulse & dec (default)
;pl5 : f2 channel - power level for 13C cpmg [6.6 kHz => 38us pulse length]
;pl6 : f4 channel - power level for 2H refocus during cpmg []
;pl9 : f1 channel - power level for 1H presat
;pl12: f2 channel - power level for 13C CPD/BB decoupling [2.5 kHz => 100us pulse length]
;sp1: f1 channel - shaped pulse H2O flipdown [default]
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse	[2*p1]
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse	[2*p3]
;p5 : f2 channel -  90 degree cpmg power pulse	[38 us]
;p6 : f2 channel - 180 degree cpmg power pulse	[2*p5]
;p9 : f1 channel - presat pulse [1-1.5s]
;p14: f1 channel - 90 degree shaped pulse [5 ms]
;p16: gradient pulse gp1, gp4 & gp5	[1 ms]
;p17: gradient pulse gp2 & gp6		[0.5 ms]
;p18: gradient pulse gp3		[0.9 ms]
;p19: gradient pulse gp7		[0.3 ms]
;p23: f4 channel -  90 degree @ pl4 [300 us]
;p25: f4 channel -  90 degree refocus @ pl6 [300 us]
;p26: f4 channel -  180 degree refocus @ pl6 [2*p25]
;d10 : incremented delay (2D)
;d1 : relaxation delay; 1-5 * T1
;d3 : total time for CPMG trains        [<= 40 ms!]
;d4 : 1/(4*JCH) taub			[1.785 ms]
;d5 : tau_cpmg
;d6 : taueq,				[400 us]
;d7 : < 1/(4*JCH) taua			[1.785 ms]
;d8 : tau_lim limiting delay 2H decoupling	[2 ms]
;d11: delay for disk I/O                [30 ms]
;d12: delay power switching             [20 us]
;d16: delay for gradient recovery       [200 us]
;d19: delay 3-9-19 watergate		[500 200 u, 600 167u, 800 125u]
;cnst1: f1 offset for presat (Hz from carrier)  [0]
;cnst3: J(CH)                           [135 Hz]
;l3: number of cpmg cycles in each period [<= 40 !]
;  total number of cycles = 2*l3
;l6: loop for 2H refocusing
;l7: l6-1, loop for 2H refocusing
;l9: flag for 2H
;spnam1: hard or Eburp2.1000
;spoffset1: 0
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2 [garp]
;cpd4: decoupling according to sequence defined by cpdprg4 [waltz16]
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence [133 us]
;pcpd4: f4 channel - 90 degree pulse for decoupling sequence [300 us]

;for z-only gradients:
;gpz1: 32.28%   10.3 G/cm       V: 5000
;gpz2: 25.82%   8.24 G/cm       V: 4000
;gpz3: 78.09%   39.1 G/cm       V: 19000
;gpz4: 71.01%   22.6 G/cm       V: 11000
;gpz5: -38.73%  12.4 G/cm       V: -6000
;gpz6: 73.98%   37.1 G/cm       V: 18000
;gpz7: 35.51%   11.3 G/cm       V: 5500

;use gradient files:
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: RECT.1
;gpnam4: SINE.100
;gpnam5: SINE.100
;gpnam6: RECT.1
;gpnam7: SINE.100

				;preprocessor-flags-start
;H2REF: for 2H refocusing pulses during cpmg. Start experiment with
;	option -DH2REF (eda: ZGOPTNS)
;H2DEC: for 2H decoupling during t1. Start experiment with
;       option -DH2DEC (eda: ZGOPTNS)
				;preprocessor-flags-end

;$Id: C13_CHD2_METmethylCPMG.gcm,v 4.0 2015/01/27 11:54:19 gcm Exp $
