;Renamed: hsqcfpf3gpphwg_t2.jk -> hsqcfpf3gpphwg_T2.jk, 3/10/2014
;
;For measuring amide proton T2s (use for PREs), Oct 2013
;Derived from hsqcfpf3gpphwg.3.2.jk
;
;J. Iwahara, C. Tang & G.M. Clore
;J. Magn. Reson. 184, 185-195 (2007)
;
;L.W. Donaldson, L.E. Kay et al
;J. Am. Chem. Soc. 123, 9843-9847 (2001)
;
;
;Option for 1D spectrum with no 15N shift evolution (-DONE_D)
;Delays adjusted for zero first-order phase correction
;Without refocusing for 0,0 phase correction
;Option for (180,-360) phase correction (default is (90,-180) phase corr.)
;Option for carbon decoupling
;Assumes p8 > p2
;With gradients during t1 to keep water along z
;With separate shape/power level for second water pulse in WATERGATE (flip-back)

;$CLASS=HighRes
;$DIM=3D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>
 
#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

; define vdlist and loop counter
define list<delay> t2delay = <$VDLIST>
"l1=0"

"p2=p1*2"
"p22=p21*2"

"d11=30m"
"d12=20u"
"d26=1s/(cnst4*4)"

"in0=inf2/4"

"DELTA1=d26-p16-d16-4u"
"DELTA2=d26-p19-d16-p11-d12-4u"
"DELTA3=d26-p19-d16-p11-d12-8u-de+0.63662*p1"

# ifndef ONE_D

# ifdef SINGLEDWELL
    "d0=in0-0.5*(10u+p1+0.63662*p21)"
# else
    "d0=in0/2-0.5*(10u+p1+0.63662*p21)"
# endif /*SINGLEDWELL*/

# endif /*ONE_D*/

# ifdef AMIDESEL
/*  REBURP (p40, sp24)   */
"cnst19=8.2"
"spoffs24=bf1*(cnst19/1000000)-o1"
"p40=2257*600.4/bf1" /* REBURP pulse length  */
"spw24=plw1*(pow((p1*1.97/p40)/0.0798,2))"   /* REBURP power level  */
"spoal24=0.5"
# endif /*AMIDESEL*/

/* 15N decoupling */
"p62=300u"
"plw26=plw3*(pow((p21/p62),2))"   /* CPD power level  */

/* 1H purge (8 kHz) */
"p6=32u"
"plw8=plw1*(pow((p1/p6),2))"   /* 8 kHz power level  */


aqseq 312


1 ze 
  d11 pl16:f3
2 d11 do:f3

# ifdef LABEL_CN
  if "4*d0+p2-24u > p8"
  {
   d12  
   "d60=d0+0.5*p1-0.25*p8"
  }
# endif /*LABEL_CN*/

  50u UNBLKGRAD
  p16:gp1		; crush all residual transverse magnetisation
  d16

# ifndef NOPROTCRUSH

;# ifndef WATERCRUSH
;    10u pl0:f1
;    (p11:sp1 ph10:r):f1	; flipdown(-x): +z -> +y
;# endif /*WATERCRUSH*/
;
;  10u pl1:f1
;  (p1 ph1):f1

; purge water before recycle delay
10u pl8:f1
2mp ph1
3mp ph2
10u pl1:f1

  4u
  p16:gp1*0.73		; crush 1H magnetisation at start of scan
  d16

# endif /*NOPROTCRUSH*/

  4u BLKGRAD
  d12
  "TAU=t2delay[l1]/4"

  d1

3 d12 pl1:f1 pl3:f3
 
  50u UNBLKGRAD
  (p21 ph1):f3
  4u
  p16:gp1*0.57
  d16*2

# ifdef AMIDESEL
  10u pl0:f1
  (p11:sp1 ph13:r):f1	; flipdown(-x/+x), +z -> +y/-y
  10u pl1:f1
# endif AMIDESEL

  (p1 ph11)		; +x/-x

;---------start relaxation/transfer period

# ifndef AMIDESEL
    5u
    5u gron0
    TAU
    5u
    5u groff
# else
    20u
    TAU
# endif /*AMIDESEL*/

  4u
  p16*0.5:gp2
  d16
  (p22 ph6):f3
  4u
  p16*0.5:gp2*-1
  d16

# ifndef AMIDESEL
    5u
    5u gron0*-1
    TAU
    5u
    5u groff
# else
    20u
    TAU
# endif /*AMIDESEL*/

  d26

  4u
  p16:gp3
  d16

# ifdef AMIDESEL
  10u pl0:f1
  (p40:sp24 ph12):f1
  10u pl1:f1
# else
  (p2 ph12)
# endif /*AMIDESEL*/ 

  4u
  p16:gp3
  d16

# ifndef AMIDESEL
    5u
    5u gron0
    TAU
    5u
    5u groff
# else
    20u
    TAU
# endif /*AMIDESEL*/

  4u
  p16*0.5:gp2
  d16
  (p22 ph6):f3
  4u
  p16*0.5:gp2*-1
  d16

# ifndef AMIDESEL
    5u
    5u gron0*-1
    TAU
    5u
    5u groff
# else
    20u
    TAU
# endif /*AMIDESEL*/

  d26

;---------end relaxation period

  (p1 ph2) 		; +y

  d12 pl0:f1
# ifdef AMIDESEL
  (p11:sp1 ph2:r):f1		; flipdown(+y), +x -> -z
# else
  (p11:sp1 ph8:r):f1		; flipdown(+x/-x), -y/+y -> -z
# endif /*AMIDESEL*/
  d12 pl1:f1

  p16:gp4
  d16

;--------------------start 15N shift evolution-------------------

# ifndef ONE_D

# ifdef LABEL_CN

  if "4*d0+p2-24u < p8"
  { 
   (p21 ph3):f3
   2u
   d0 gron0
   d0 gron0*-1
   8u groff
   (p2 ph5):f1
   2u
   d0 gron0
   d0 gron0*-1
   8u groff
   (p21 ph4):f3
  }
  else
  {
  (p21 ph3):f3
   2u
   d60 gron0
   d60 gron0*-1
   8u groff
   ( center (p2 ph5):f1 (p8:sp13 ph1):f2 )
   2u
   d60 gron0
   d60 gron0*-1
   8u groff
   (p21 ph4):f3
   }

# else

  (p21 ph3):f3
  2u
  d0 gron0
  d0 gron0*-1
  8u groff
  (p2 ph5):f1
  2u
  d0 gron0
  d0 gron0*-1
  8u groff
  (p21 ph4):f3

# endif /*LABEL_CN*/

# else

  (center (p2 ph5):f1 (p21 ph3 6u p21 ph4):f3 )

# endif /*ONE_D*/

;---------------------end 15N shift evolution--------------------

  4u
  p16:gp6
  d16 pl0:f1
  (p11:sp1 ph7:r):f1	; flipdown(-x), +z -> +y
  4u
  4u pl1:f1

  (p1 ph1)
  4u
  p19:gp5
  d16 pl0:f1
  DELTA2
  (p11:sp1 ph7:r):f1    ; flipdown(-x), +z -> +y
  d12 pl1:f1
  (center (p2 ph1) (p22 ph1):f3 )
  d12 pl0:f1
  (p11:sp11 ph9:r):f1   ; flipback(-x), -y -> +z
  4u
  p19:gp5
  d16 pl26:f3
  DELTA3
  4u BLKGRAD

  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
	F1QF(iu1)
	F2PH(ip3, id0)
exit 
  

ph1= 0
ph2= 1
ph3= 0 2
ph4= 0 0 0 0 2 2 2 2
ph5= 0 0 0 0 2 2 2 2
ph6= 0
ph7= 2
ph8= 0 0 2 2
ph9= 2
ph10=2
ph11=0 0 2 2
ph12=1 1 3 3
ph13=2 2 0 0
ph31=0 2 2 0 2 0 0 2


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl8: 1H purge power (ca. 10 kHz) 
;pl26: f3 channel - power level for CPD/BB decoupling
;sp1: f1 channel - shaped pulse  90 degree (sinc flip-down)
;sp11: f1 channel - shaped pulse  90 degree (sinc flip-back)
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;sp24: f1 channel - shaped pulse 180 degree (Reburp.1000)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p11: f1 channel -  90 degree shaped pulse [1000 usec]
;p16: homospoil/gradient pulse
;p19: second homospoil/gradient pulse
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p40: f1 channel - 180 degree shaped pulse for refocusing
;                      Reburp.1000       (1679/2257 usec at 600 MHz)
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d26 : 1/(4J)YH
;cnst4: = J(YH)
;cnst19: H(N) chemical shift (offset, in ppm)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd3: decoupling according to sequence defined by cpdprg3
;cpdprg3: garp4.p62
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz0: 1-2%
;gpz1: 73.73%
;gpz2: 7%
;gpz3: 11%
;gpz4: 17%
;gpz5: 53%
;gpz6: 13%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end

