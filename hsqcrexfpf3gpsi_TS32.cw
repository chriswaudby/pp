;don't use this one - doesn't seem to work!
;
;15N constant-time HSQC CPMG relaxation dispersion experiment
;with 1H CW decoupling during CPMG period, pseudo-3D (recommended version)
;Only for use with Topspin 3
;Based on hsqcrexfpf3gpsi.2.jk 
;
;D.F. Hansen, P. Vallurupalli & L.E. Kay
;J. Phys. Chem. B 112, 5898-5904 (2008)
;
;John Kirkpatrick, UCL, Feb 2014
;
;See end of sequence for available options
;Using external power list for proton CW
;Proton CW compensation for reference plane, 15N compensation for all planes
;15N compensation block after recycle delay and immediately before main sequence
;1H compensation block for reference plane after FID
;Nb. DELTA5 can become negative for long 1H pw90 and high 1H CW field-strength
;
;------------------------------------------------------------------------------------------
;
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<power> protCW= <$VALIST>		; list of power levels for proton CW

define loopcounter COUNTER1

"p2=p1*2"
"p22=p21*2"
"p25=90u"		; 180deg pulse for CPMG

"TAU3=0.6366*p21"	; delay equal to phase evolution during 90deg CPMG-flanking pulses
"TAU4=0.3634*p21"

# ifdef CAL_N
    define pulse pNcal
    "pNcal=p25/2"
# endif /*CAL_N*/

"plw23=plw3*(p22/p25)*(p22/p25)*cnst63"		; power for CPMG train

"d11=30m"
"d12=20u"
"d13=4u"

"d0=10u"
"in0=inf2/4"

"d26=1s/(cnst4*4)"

# ifdef LABEL_CN
    "DELTA=4*d0+larger(p8,2*p2)+26u"
# else
    "DELTA=4*d0+2*p2+26u"
# endif /*LABEL_CN*/

"DELTA1=d25-p16-d16-4u-larger(p1,p21)"
"DELTA2=d26-p19-d16-4u-larger(p1,p21)"
"DELTA3=d26-p16-d16-4u-larger(p1,p21)"
"DELTA4=p16+d16+8u+de-0.6366*p1"

"DELTA6=d12+0.3634*p21-0.6366*p1"

aqseq 312

1 ze
  d11

# ifndef HPDISP

  if "d1 < 1.0s"
    {
    d11
    print "error: recycle delay too short, d1 must be >= 2.0 s, aborting..."
    goto stop
    }

# else /*HPDISP*/

  d11
  print "warning: hpdisp option is turned on; turn off for genuine acquisition (and check d1)."

# endif /*HPDISP*/

  if "d21 > 25m"
    {
    d11
    print "error: CPMG train too long, d21 must be =< 25 ms, aborting..."
    goto stop
    }

  if "p25 < 80u"
    {
    d11
    print "error: CPMG 180deg pulse too short, p25 must be >= 80 us, aborting..."
    goto stop
    }


# if defined (INTERLEAVE) && !defined (HPDISP)
    d11 st0
# else /*INTERLEAVE, HPDISP*/
    define loopcounter PHASE
    define loopcounter VD_CNT
    define loopcounter T1_CNT
    d12
    "VD_CNT=td1"
    if "td2 == 1"
      {
      d12
      "PHASE=1; T1_CNT=1"
      }
    else
      {
      d12
      "PHASE=2; T1_CNT=td2/2"
      }
# endif /*INTERLEAVE, HPDISP*/

# ifdef LABEL_CN
    d12 pl0:f2
# endif /*LABEL_CN*/

  d12
;  "protCW.idx=0"	; this statement should work, but doesn't!
  3m protCW.res

2 30m
3 d11
4 9m
5 d13 do:f3

  d12
  "l10=vd"

  if "l10 > 1500"
    { 
    d11
    print "error: CPMG frequency too high, must be < 1500, aborting..."
    goto stop
    }

  if "l10 == 0"
    {
    d12*2
    "TAU=1m; TAU1=1m; COUNTER=0"		; precaution
    d12*2
    "p32=250000/cnst52; d32=p32"		; equivalent pw90 for 1H CW
    d12
    }
  else
    {
    d12
    "TAU=(1s / (l10*4) ) - p25/2"		; NOTE SYNTAX
    d12
    "TAU1=TAU-TAU3-10u"
    d12
    "COUNTER=d21*l10*2 - 0.9"			; NOTE MINUS SIGN
    d12
    "FACTOR1=cnst52/(2*l10)"			; multiplier for 1H CW field-strength
    d12
    "p32=250000/(FACTOR1*l10*2); d32=p32"	; equivalent pw90 for 1H CW
    }

  d12
  "DELTA5=0.6366*p32-1.27324*p1"

  d1	;---------------------------------------- recycle delay

# ifndef NOCOMP

  d12
  "l11=cnst62-vd"				; for 15N compensation block

  if "cnst62 > 1500"
    { 
    d11
    print "error: maximum CPMG frequency as specified by cnst62 is too high, must be < 1500, aborting..."
    goto stop
    }

  if "l11 < 0"
    {
    d11
    print "error: cnst62 is less than maximum CPMG frequency in vdlist, aborting..."
    goto stop
    }

  if "l11 == 0"
    {
    d12*2
    "TAU2=1s; COUNTER1=0"			; precaution
    d12
    d21*2
    }
  else
    {
    d12
    "TAU2=(1s / (l11*4) ) - p25/2"		; NOTE SYNTAX
    d12
    "COUNTER1=d21*l11*2 + 0.1"			; NOTE PLUS SIGN
    d12 pl23:f3
6   TAU2
    (p25 ph2):f3
    TAU2
    TAU2
    (p25 ph2):f3
    TAU2
    lo to 6 times COUNTER1
    }

# endif /*NOCOMP*/


  50u UNBLKGRAD
  d12 pl3:f3
  p16:gp6*0.7
  d16
  (p21 ph1):f3
  4u
  p16:gp6
  d16 

  d12 pl0:f1
  (p11:sp1 ph4:r):f1		; flipdown(-y), +z -> -x
  10u pl1:f1
  
  (p1 ph1):f1
  4u
  p16:gp1
  d16
  DELTA1
  (center (p2 ph1):f1 (p22 ph1):f3 )
  DELTA1
  p16:gp1
  d16
  (p1 ph2):f1
  4u
  p16:gp2
  d16

  (p21 ph1):f3
  4u
  p19:gp3
  d16
  DELTA2
  (center (p2 ph1):f1 (p22 ph1):f3 )
  4u
  p19:gp3
  d16
  DELTA2
  p21:f3 ph2 p2:f1 ph1
  4u
  p16:gp4
  d16

  d8				; re-equilibration delay

;-----------------start CPMG

  30u fq=cnst19(bf ppm):f1	; move proton offset to amides
  
  (p1 ph11):f1
  DELTA5
  (p1 ph1):f1
  DELTA6
  (p2 ph11):f1
  d12 protCW:f1
  
  if "l10 == 0"
    {
    (p21 ph12):f3
    d13 10u
    }
  else
    {
    (ralign (TAU3 cw ph1):f1 (p21 ph12):f3 )
    TAU1
    10u pl23:f3
    (p25 ph2):f3
    TAU
7   TAU
    (p25 ph2):f3
    TAU
    lo to 7 times COUNTER
    d13 do:f1
    10u pl3:f3
    }

  10u pl1:f1
  (center (p2 ph1):f1 (p22 ph13):f3 )
  10u protCW:f1

  if "l10 == 0"
    {
    d13 10u
    (p21 ph1):f3
    }
  else
    {
    d13
    10u pl23:f3
8   TAU cw:f1 ph1
    (p25 ph2):f3
    TAU
    lo to 8 times COUNTER
    TAU
    (p25 ph2):f3
    10u pl3:f3
    TAU1
    (lalign (TAU3 TAU4 do):f1 (p21 ph1):f3 )
    }

  d12 pl1:f1
  (p2 ph14):f1
  DELTA6
  (p1 ph1):f1
  DELTA5
  (p1 ph14):f1
  
  30u fq=0:f1			; move proton offset back to water

;-----------------end CPMG

# ifdef CAL_N			; optional calibration element
    p16:gp12*0.7
    d16 pl23:f3
    (pNcal ph1):f3
    4u
    p16:gp12
    d16 pl3:f3
# endif /*CAL_N*/

  d8				; re-equilibration delay

  p16:gp5
  d16

;-----------------start 15N chemical shift evolution

  (p21 ph5):f3
  4u
  d0 gron0
  d0 gron0*-1
  10u groff 
  p16:gp7*EA
  d16

# ifdef LABEL_CN
    (center (p1 ph21 1u p2 ph22 1u p1 ph21):f1 (p8:sp13 ph1):f2 )
# else
    (p1 ph21 1u p2 ph22 1u p1 ph21):f1	; composite 1H inversion
# endif /*LABEL_CN*/

  4u
  d0 gron0
  d0 gron0*-1
  10u groff
  (p22 ph1):f3
  4u
  p16:gp7*EA*-1
  d16
  DELTA

;-----------------end 15N chemical shift evolution

  (center (p1 ph3):f1 (p21 ph6):f3 )
  4u
  p16:gp8
  d16
  DELTA3
  (center (p2 ph1):f1 (p22 ph1):f3 )
  4u
  DELTA3
  p16:gp8
  d16
  (center (p1 ph2):f1 (p21 ph2):f3 )
  4u
  p16:gp9
  d16
  DELTA3
  (center (p2 ph1):f1 (p22 ph1):f3 )
  4u
  DELTA3
  p16:gp9
  d16
  (p1 ph1):f1
  DELTA4
  (p2 ph1):f1
  4u
  p16:gp10
  d16 pl16:f3
  4u

# ifdef INTERLEAVE
    goscnp ph31 cpd3:f3
# else
    gosc ph31 cpd3:f3
# endif /*INTERLEAVE*/

  4u do:f3

# ifndef NOCOMP

  if "l10 == 0"
    {
    d12 protCW:f1
    (p32 ph2):f1
    d13
    d13 cw:f1 ph1
    d21*2
    d13 do:f1
    (p32 ph4):f1
    }
  else
    {
    d12 d32*2 d13*3
    }

  p16:gp11*0.7
  d16 pl0:f1
  (p11:sp1 ph3:r):f1		; flipdown(-x), +z -> +y
  10u pl1:f1
  (p1 ph1):f1
  p16:gp11
  d16

# endif /*NOCOMP*/

  4u BLKGRAD

# ifdef INTERLEAVE
  
# ifndef HPDISP

  3m st ivd
  3m protCW.inc
  lo to 2 times nbl
  15m protCW.res
  15m ipp5 ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp31
  lo to 3 times ns

  d11 mc #0 to 3
     F1QF()
     F2EA(calgrad(EA) & calph(ph6, +180) & exec(rppall), caldel(d0, +in0) & calph(ph5, +180) & calph(ph31, +180))

# else /*HPDISP*/

  3m ivd
  3m protCW.inc
  lo to 2 times VD_CNT
  15m protCW.res
  15m ipp5 ipp11 ipp12 ipp13 ipp14 ipp21 ipp22 ipp31
  lo to 3 times ns

  10m wr #0 if #0 zd
  10m igrad EA 5m ip6*2 5m rpp5 rpp11 rpp12 rpp13 rpp14 rpp21 rpp22 rpp31
  lo to 4 times PHASE
  3m id0 3m ip5*2 3m ip31*2
  lo to 5 times T1_CNT

# endif /*HPDISP*/

# else /*INTERLEAVE*/

  6m
  lo to 2 times ns
  10m wr #0 if #0 zd 10m ivd 10m protCW.inc
  lo to 3 times VD_CNT
  10m protCW.res 10m igrad EA 10m ip6*2
  lo to 4 times PHASE
  3m id0 3m ip5*2 3m ip31*2
  lo to 5 times T1_CNT
  
;;------------------------ this construct is not working in Topspin 3.2.5 --------------------------
;  d11 mc #0 to 3 
;     F1QF(ivd & protCW.inc)
;     F2EA(exec(protCW.res) & calgrad(EA) & calph(ph6, +180), caldel(d0, +in0) & calph(ph5, +180) & calph(ph31, +180))
;;--------------------------------------------------------------------------------------------------

# endif /*INTERLEAVE*/

stop, exit


ph1= 0 
ph2= 3
ph3= 2
ph4= 1
ph5= 3 3 3 3 1 1 1 1
ph6= 0
ph11=3 3 1 1
ph12=0 2
ph13=0 0 2 2
ph14=1 1 3 3
ph21=0 0 2 2
ph22=3 3 1 1
ph31=0 2 0 2 2 0 2 0
   
;----Varian phases

;ph1= 0 
;ph2= 1
;ph3= 2
;ph4= 3
;ph5= 1 1 1 1 3 3 3 3	; phi5
;ph6= 0			; phi6
;ph11=1 1 3 3 		; phi1
;ph12=0 2		; phi2
;ph13=0 0 2 2		; phi3
;ph14=3 3 1 1		; phi4
;ph21=0 0 2 2
;ph22=1 1 3 3
;ph31=0 2 0 2 2 0 2 0
  

;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl16: f3 channel - power level for CPD/BB decoupling
;pl23: f3 channel - power level for CPMG (calculated internally)
;sp1:  f1 channel - shaped pulse  90 degree [flipdown]
;sp11: f1 channel - shaped pulse  90 degree [flipback]
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse                         [1 msec]
;p19: gradient pulse 2                                 [2 msec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p25: f3 channel - 180 degree pulse at pl23            [90 usec]
;p32: f1 channel -  90 degree pulse for proton CW (calculated internally)
;d0 : incremented delay for 15N chemical shift evolution
;d1 : recycle delay; 1-5 * T1
;d8 : equilibration delay                              [2-5 msec]
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d13: short delay                                      [4 usec]
;d16: delay for homospoil/gradient recovery
;d21: length of single CPMG train
;     total constant-time relaxation period is 2*d21
;d25: =< 1/(4J(NH))
;d26: 1/(4J(NH))
;d32: delay equivalent to pulse length p32 (set internally)
;vd : variable field strength, taken from vd-list
;cnst4: = J(NH)
;cnst19: offset for centre of amides in ppm
;cnst50: power level for proton hard pulse (set internally)
;cnst51: power level for proton CW (calculated internally)
;cnst52: nominal field strength for proton CW in Hz
;cnst62: maximum CPMG frequency in Hz
;cnst63: additive (multiplicative) correction to power-level (power) for 15N CPMG (calibrate)
;l10: current CPMG frequency (calculated internally)
;l11: diff. between max. and current CPMG frequency (calculated internally)
;inf1: 1/SW(N) = 4 * DW(N)
;in0: 1/(4 * SW(N)) = DW(N)
;nd0: 4
;NS: 8 * n
;DS: >= 8
;td1: number of delays in vd-list
;td2: number of experiments in F2
;NBL: = td1
;FnMODE: QF in F1
;FnMODE: echo-antiecho in F2
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz0:  1-2%
;gpz1:  7%
;gpz2: -17%
;gpz3:  11.5%
;gpz4:  13%
;gpz5:  37%
;gpz6:  19%
;gpz7:  80%
;gpz8:  5%
;gpz9:  3%
;gpz10: 16.2% 
;gpz11: 23% 
;gpz12: 53% 

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam7: SMSQ10.100
;gpnam8: SMSQ10.100
;gpnam9: SMSQ10.100
;gpnam10: SMSQ10.100
;gpnam11: SMSQ10.100
;gpnam12: SMSQ10.100

;note: the values in the vd-list are interpreted as field-strength in Hz

                                          ;preprocessor-flags-start
;INTERLEAVE: for looping over vd-list before ns
;LABEL_CN: for 13C decoupling during indirect 15N evolution period
;CAL_N: for calibration of 180deg pulse for 15N CPMG
;NOCOMP: for running without compensation elements
;HPDISP: for pulse sequence display with hpdisp tool
                                          ;preprocessor-flags-end

