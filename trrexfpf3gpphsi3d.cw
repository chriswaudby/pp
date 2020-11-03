;15N constant-time TROSY CPMG relaxation dispersion with temperature compensation, pseudo-3D (recommended version)
;Only for use with Topspin3
;Based on trrexfpf3gpphsi3d.jk
;
;P. Vallurupalli, D.F. Hansen, E. Stollar, E. Meirovitch & L.E. Kay,
;PNAS 104, 18473-18477 (2007)
;
;John Kirkpatrick, UCL, Feb 2014
;
;Temperature compensation element after recycle delay
;See end of sequence for available options
;Added element to check that d21 and field strengths in vdlist file are all compatible, 12/5/14
;Removed assignment of vd entries to loop counters (to allow non-integer field strengths), 14/5/14
;
;-----------------------------------------------------------------------------------
;
;
;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple2>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p22=p21*2"

# ifdef CAL_N
    define pulse pNcal
    "pNcal=p25/2"
# endif /*CAL_N*/

"d11=30m"
"d12=20u"
"d13=4u"

"d26=1s/(cnst4*4)"

# ifndef ONE_D
"in10=inf2/4"
"d10=in10/2-5u"
# endif /*ONE_D*/

define delay VDCOMP
define loopcounter COUNTER1

;----parameters needed to check that CPMG counter is always naturally integer for given d21/vdlist

define loopcounter VDMIN
define loopcounter D21_USEC

"D21_USEC=d21*1000000 + 0.5"	; length of d21 in microseconds (integer)
"VDMIN=500000000/D21_USEC"	; minimum allowed CPMG field-strength x 1000 (integer)

;----

"DELTA1=d25-p16-d16-larger(p1,p21)-4u"
"DELTA2=d26-p16-d16+2u"
"DELTA3=d25-p19-d16-larger(p1,p21)-d12-p11-p1"
"DELTA4=d25-p19-d16-larger(p1,p21)-d12-p11-4u"
"DELTA5=d25-p19-d16-larger(p1,p21)-d12-p11-p21-4u"

"plw23=plw3*(p22/p25)*(p22/p25)*cnst63"		; power for CPMG train

# ifndef ANTI_TROSY
    "TAU2=4*p1+4u"
# endif /*ANTI_TROSY*/

"l0=1"

"spoff1=0"
"spoff11=0"

# ifdef LABEL_CN
    "spoff13=0"
# endif /*LABEL_CN*/


aqseq 312


1 ze
  d11

  if "d1 < 2.0s"
    {
    d11
    print "error: recycle delay too short, d1 must be >= 2.0 s, aborting..."
    goto stop
    }

  if "d21 > 25.1m"
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

# ifdef INTERLEAVE
    d11 st0
# else
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
# endif /*INTERLEAVE*/

2 20m
3 d11
4 12m
5 d13

# if defined(LABEL_CN) && !defined (ONE_D)
    if "4*d10-24u > p8"
      {
      d12
      "d60=d10-0.25*p8"
      }
    d12 pl0:f2
# endif /*LABEL_CN, ONE_D*/


  if "vd > 1500"
    {
    d11
    print "error: CPMG frequency too high, must be < 1500, aborting..."
    goto stop
    }

  if "(1000*vd)%VDMIN != 0"
    {
    d11
    print "error: One or more field strengths in vdlist file does not give integer loop counter \
           for CPMG train - check vdlist and/or d21.  Aborting..."
    goto stop
    }

  if "vd == 0"
    {
    20u
    "TAU=1s"
    20u
    "COUNTER=0"
    }
  else
    {
    20u
    "TAU=(1 / (vd*4) ) - p25/2000000"
    20u
    "COUNTER=d21*vd*2 + 0.1"
    }

  d1

# ifndef NOCOMP

  20u
  "VDCOMP=cnst62-vd"

  if "cnst62 > 1500"
    {
    d11
    print "error: maximum CPMG frequency as specified by cnst62 is too high, must be < 1500, aborting..."
    goto stop
    }

  if "VDCOMP < 0"
    {
    d11
    print "error: cnst62 is less than maximum CPMG frequency in vdlist, aborting..."
    goto stop
    }

  if "VDCOMP == 0"
    {
    20u
    "TAU1=1s"
    20u
    "COUNTER1=0"
    }
  else
    {
    20u
    "TAU1=(1 / (VDCOMP*4) ) - p25/2000000"
    20u
    "COUNTER1=d21*VDCOMP*2 + 0.1"
    }

  if "VDCOMP == 0"
    {
    d12
    d21*2
    }
  else
    {
    if "vd == 0"
      {
      d21*2
      }
    d12 pl23:f3
6   TAU1
    (p25 ph21):f3
    TAU1
    TAU1
    (p25 ph21):f3
    TAU1
    lo to 6 times COUNTER1
    }

# endif /*NOCOMP*/

  50u UNBLKGRAD
  d12 pl3:f3
  p16:gp0*0.7
  d16
  (p21 ph20):f3
  4u
  p16:gp0
  d16*2

  d12 pl0:f1
  (p11:sp1 ph13:r):f1		; flipdown(+y): +z -> +x
  d12 pl1:f1

# ifdef ANTI_TROSY
    (p1 ph10)
# else
    (p1 ph12)
# endif /*ANTI_TROSY*/

  4u
  p16:gp1
  d16
  DELTA1
  (center (p2 ph10) (p22 ph20):f3 )
  DELTA1
  4u
  p16:gp1
  d16
  (p1 ph11)
  
  4u
  p16:gp2
  d16*2
  d12 pl23:f3

  d8			; tau_eq

  (p25*0.5 ph1):f3

  if "vd > 0"
    {
7   TAU
    (p25 ph21):f3
    TAU
    lo to 7 times COUNTER
    }

  (p25*0.5 ph2):f3

  DELTA2
  p16:gp3
  d16 pl3:f3 

  (p1 ph11 2u p2 ph10 2u p1 ph11):f1
  2u  

# ifdef ANTI_TROSY
    (p21 ph20 2u p22 ph21 2u p21 ph20):f3	; 180(y)
# else
    (p21 ph21 2u p22 ph20 2u p21 ph21):f3	; 180(x)
# endif /*ANTI_TROSY*/

  DELTA2 
  p16:gp3
  d16 pl23:f3

  (p1 ph11 2u p2 ph10 2u p1 ph11):f1 
  2u
  (p25*0.5 ph3):f3

  if "vd > 0"
    {
8   TAU
    (p25 ph20):f3
    TAU
    lo to 8 times COUNTER
    }

  (p25*0.5 ph21):f3

# ifdef CAL_N
    4u
    p16:gp8*0.7
    d16
    (pNcal ph20):f3
    4u
    p16:gp8
    d16
# endif /*CAL_N*/

  d8			; tau_eq

# ifdef ANTI_TROSY
    (p1 ph11 2u p2 ph10 2u p1 ph11):f1
# else
    TAU2
# endif /*ANTI_TROSY*/

  4u
  p16:gp4
  d16*2
  d12 pl3:f3

  if "l0 %2 == 1"
    {
    (p21 ph4):f3
    }
  else
    {
    (p21 ph14):f3
    }

# ifndef ONE_D

# ifdef LABEL_CN

    if "4*d10-24u < p8"
      {
      2u
      d10 gron5
      d10 gron5*-1
      10u groff
      d10 gron5
      d10 gron5*-1
      8u groff
      }
      else
      {
      2u
      d60 gron5
      d60 gron5*-1
      8u groff
      (p8:sp13 ph10):f2
      2u
      d60 gron5
      d60 gron5*-1
      8u groff
      }

# else /*LABEL_CN*/
  
    2u
    d10 gron5
    d10 gron5*-1
    10u groff
    d10 gron5
    d10 gron5*-1
    8u groff

# endif /*LABEL_CN*/

# endif /*ONE_D*/

  (p1 ph5)

  d12 pl0:f1

# ifdef ANTI_TROSY
    (p11:sp11 ph5:r):f1		; flipback(+y): -x -> +z
# else
    (p11:sp11 ph6:r):f1         ; flipback(-y): +x -> +z
# endif /*ANTI_TROSY*/

  DELTA3
  p19:gp6
  d16 pl1:f1

  (center (p2 ph10) (p22 ph20):f3 )

  4u
  p19:gp6
  d16 pl0:f1 
  DELTA4

  (p11:sp11 ph10:r):f1		; flipback(+x): -z -> +y
  d12 pl1:f1

  (center (p1 ph10) (p21 ph7):f3 )

  4u
  p19:gp7
  d16 pl0:f1
  DELTA4

  (p11:sp1 ph12:r):f1		; flipdown(-x): +z -> +y
  d12 pl1:f1

  (center (p2 ph10) (p22 ph20):f3 )

  d12 pl0:f1
  (p11:sp11 ph22:r):f1		; flipback(-x): -y -> +z

  DELTA5
  p19:gp7
  d16 pl1:f1 

  4u BLKGRAD
  (p21 ph20):f3

# ifdef INTERLEAVE

  goscnp ph31
  3m st ivd
  lo to 2 times nbl

  20m ipp1 ipp2 ipp3 ipp4 ipp14 ipp31
  lo to 3 times ns

  d11 mc #0 to 3
     F1QF()
     F2EA(calph(ph5, +180) & calph(ph6, +180) & calph(ph7, +180) & calph(ph31, +180) & calclc(l0, 1) & exec(rppall), caldel(d10, +in10) & calph(ph4, +180) & calph(ph14, +180) & calph(ph31, +180))

# else /*INTERLEAVE*/

  gosc ph31
  3m
  lo to 2 times ns
  10m wr #0 if #0 zd 10m ivd
  lo to 3 times VD_CNT
  6m ip5*2 6m ip6*2 6m ip7*2 6m ip31*2 6m iu0
  lo to 4 times PHASE
  3m id10 3m ip4*2 3m ip14*2 3m ip31*2
  lo to 5 times T1_CNT

;;--------------------------- this construct is not working in Topspin 3.2.5  ---------------------------
;  d11 mc #0 to 3
;     F1QF(ivd)
;     F2EA(calph(ph5, +180) & calph(ph6, +180) & calph(ph7, +180) & calph(ph31, +180) & calclc(l0, 1), caldel(d10, +in10) & calph(ph4, +180) & calph(ph14, +180) & calph(ph31, +180))
;;-------------------------------------------------------------------------------------------------------

# endif /*INTERLEAVE*/

stop, exit
   

; Proton/carbon phases

ph5= 1
ph6= 3
ph10=0
ph11=3
ph12=2
ph22=2
ph13=1

; Nitrogen phases

ph1= 0 2
ph2= 3 3 1 1 
ph3= 0 0 2 2
ph4= 3 3 1 1 2 2 0 0
ph14=3 3 1 1 0 0 2 2
ph7= 1
ph20=0
ph21=3

; Receiver phase

ph31=3 1 1 3 0 2 2 0


;;-------------Varian phases----------------
;
;; Proton phases
;
;ph5= 3
;ph6= 1
;ph10=0
;ph11=1
;ph12=2
;ph13=3
;
;; Nitrogen phases
;
;ph1= 0 2
;ph2= 1 1 3 3 
;ph3= 0 0 2 2
;ph4= 1 1 3 3 2 2 0 0
;ph14=1 1 3 3 0 0 2 2
;ph7= 3
;ph20=0
;ph21=1
;
;; Receiver phase
;
;ph31=1 3 3 1 0 2 2 0
;
;;------------------------------------------
  

;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl23: f3 channel - power level for CPMG
;sp1:  f1 channel - shaped pulse  90 degree [flipdown]
;sp11: f1 channel - shaped pulse  90 degree [flipback]
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse                       [1 msec]
;p19: gradient pulse 2                               [500 usec]
;p21: f3 channel -  90 degree high power pulse
;pNcal: f3 channel - calibration 90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p25: f3 channel - 180 degree pulse at pl23          [90 usec]
;d1 : relaxation delay; 1-5 * T1
;d8 : equilibration delay
;d10: incremented delay                              [3 usec]
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21: length of single CPMG train
;     total constant-time relaxation period is 2*d21
;d25: 1/(4J)YH for YH
;d26: 1/(4J(YH))
;d60: calculated internally
;vd : variable field strength, taken from vd-list
;cnst4: = J(YH)
;cnst62: maximum CPMG frequency in vd-list
;cnst63: multiplicative correction to power for 15N CPMG (calibrate)
;vd: taken from vdlist (CPMG field strength)
;VDCOMP: calculated internally (CPMG field strength for compensation element)
;inf2: 1/SW(X) = 4 * DW(X)
;in10: 1/(4 * SW(X)) = DW(X)
;nd10: 4
;NS: 8 * n
;DS: >= 32
;td1: number of delays in vd-list
;td2: number of experiments in F2
;NBL: = td1
;FnMODE: QF in F1
;FnMODE: echo-antiecho in F2


;for z-only gradients:
;gpz0: -37%
;gpz1: 11%
;gpz2: 29%
;gpz3: 7%
;gpz4: 19%
;gpz5: 1-2%
;gpz6: 23%
;gpz7: 43%
;gpz8: 71%

;use gradient files:   
;gpnam0: SINE.100
;gpnam1: SINE.100
;gpnam2: SINE.100
;gpnam3: SINE.100
;gpnam4: SINE.100
;gpnam6: SINE.100
;gpnam7: SINE.100
;gpnam8: SINE.100


;note: the values in the vd-list are interpreted as field-strength in Hz

                                          ;preprocessor-flags-start
;INTERLEAVE: for looping over vd-list before ns
;LABEL_CN: for 13C decoupling during indirect 15N evolution period
;CAL_N: for calibration of 180deg pulse for 15N CPMG
;ONE_D: for running 1D or pseudo-2D version (no 15N chemical shift evolution)
;NOCOMP: for running without compensation elements
;ANTI_TROSY: for measuring anti-TROSY line during CPMG
                                          ;preprocessor-flags-end

