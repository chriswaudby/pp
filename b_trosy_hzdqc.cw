;BEST-TROSY-H(Z/D)QC
;Chris Waudby, June 2018
;
;options:
; -DDQ = HDQC (otherwise runs HZQC)
; -DONE_D = first-row
; -DOFFRES_PRESAT = presat, pl9 on cnst21 (Hz bf)

prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"d11=30m"
"d12=20u"
"d13=4u"
"d21=1s/(cnst4*4)"

"p22=p21*2"

"in0=inf1"
# ifdef ONE_D
"d0=2u"
#else
"d0=in0/2-p21*4/3.1415"
# endif /*ONE_D*/

"d2=p39-p39*cnst39-0.3633*p21"
"d3=0.5*p40-0.3633*p21"
"DELTA1=d21-p39*cnst39-p40*0.5-p16-d16-4u"
"DELTA2=d21-0.3633*p21-p16-d16-4u-0.5*p40"
"DELTA3=d21-p40-p16-d16-4u"
"DELTA4=d21-0.5*p40-p16-d16-4u-p21-de"
"acqt0=de"

#ifdef LABEL_CN
"d10=DELTA3+d3+p21+d0*0.5-p8*0.5"
"d9=DELTA3+d3+p21+d0*0.5-p8*0.5"
"in10=in0*0.5"
"in9=in0*0.5"
#endif

# ifdef OFFRES_PRESAT
  "TAU=d1-d11-60u-d12*2-d13-d12-50u-p21-2*p16-2*d16-12u"
# else
  "TAU=d1-d11-d12-50u-p21-2*p16-2*d16-12u"
# endif /*OFFRES_PRESAT*/


;"spoff23=bf1*(cnst19/1000000)-o1"
;"spoff24=bf1*(cnst19/1000000)-o1"
"spoff23=0" ; for amides on-resonance (recommended)
"spoff24=0"


"l0=1"  ; loop counter for shifting 1H 180 pulse between echo/anti-echoes


1 ze 
  d11 
2 d11 

  4u UNBLKGRAD
  p16:gp3
  d16
  4u BLKGRAD

# ifdef OFFRES_PRESAT
  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  TAU cw:f1 ph1
  d13 do:f1
  d12 pl1:f1
  30u fq=0:f1
# else
  TAU
# endif /*OFFRES_PRESAT*/

3 d12 pl3:f3
  50u UNBLKGRAD

  ; purge Nz
  (p21 ph1):f3
  4u
  p16:gp0
  d16

  ; begin main sequence
  if "l0 %2 == 1"
     {
  ;(lalign (p39:sp23 ph10) (d2 p21 ph11):f3 ) 
  (p39:sp23 ph10) (d2 p21 ph11):f3
     }
  else
     {
  ;(lalign (p39:sp23 ph10) (d2 p21 ph21):f3 ) 
  (p39:sp23 ph10) (d2 p21 ph21):f3
     }

  DELTA1
  4u
  p16:gp1
  d16
  (center (p40:sp24 ph1) (p22 ph12):f3 )
  4u
  p16:gp1
  d16


  if "l0 %2 == 1"
     {
#ifdef LABEL_CN
  (ralign (p40:sp24 ph16 DELTA3) (DELTA2 p21 ph13 d0 p21 ph1 d3 DELTA3):f3 (p8:sp13 ph1 d10):f2 )
#else
  (ralign (p40:sp24 ph16) (DELTA2 p21 ph13 d0 p21 ph1 d3):f3 )
  DELTA3
#endif /*LABEL_CN*/
     }
  else
     {
#ifdef LABEL_CN
  (DELTA3 p40:sp24 ph16) (DELTA3 d3 p21 ph23 d0 p21 ph1 DELTA2):f3 (d9 p8:sp13 ph1):f2
#else
  DELTA3
  (p40:sp24 ph16) (d3 p21 ph23 d0 p21 ph1 DELTA2):f3 
#endif /*LABEL_CN*/
     }
  4u
  p16:gp2
  d16
  (center (p40:sp24 ph1) (p22 ph1):f3 )
  4u
  p16:gp2
  d16
  DELTA4 BLKGRAD
  (p21 ph14):f3

  go=2 ph31 
#ifdef LABEL_CN
  d11 mc #0 to 2 
     F1EA(iu0 & ip13*2 & ip14*2, id0 & id10 & id9 & ip10*2 & ip31*2)
#else
  d11 mc #0 to 2 
     F1EA(iu0 & ip13*2 & ip14*2, id0 & ip10*2 & ip31*2)
#endif /*LABEL_CN*/

exit 
  
ph1=0 
ph10=0
ph11=2 0 3 1 
ph21=2 0 1 3
ph12=0 
#ifdef DQ
ph13=1 3 2 0
ph23=1 3 0 2
ph14=3
#else /* ZQ */
ph13=1 3 0 2
ph23=1 3 2 0
ph14=1
#endif
ph16=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph31=0 2 3 1 2 0 1 3

;pl3 : f3 channel - power level for pulse (default)
;pl26: f3 channel - power level for CPD/BB decoupling (low power)
;sp23: f1 channel - shaped pulse 120 degree 
;                   (Pc9_4_120.1000 or Q5.1000)
;sp24: f1 channel - shaped pulse 180 degree (Rsnob.1000)
;p16: homospoil/gradient pulse                       [1 msec]
;p21: f3 channel -  90 degree high power pulse
;p39: f1 channel - 120 degree shaped pulse for excitation
;                      Pc9_4_120.1000 (120o)    (3.0ms at 600.13 MHz)
;                  (or Q5.1000 (90o)            (2.0ms at 600.13 MHz) )
;p40: f1 channel - 180 degree shaped pulse for refocussing
;                      Rsnob.1000               (1.0ms at 600.13 MHz)
;d0 : incremented delay (2D) = in0/2-p21*4/3.1415
;d1 : relaxation delay
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d21 : 1/(4J)NH
;cnst4: = J(NH)
;cnst19: H(N) chemical shift (offset, in ppm)
;cnst39: compensation of chemical shift evolution during p39
;           Pc9_4_90.1000: 0.514
;           Pc9_4_120.1000: 0.529
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/ SW(N) = 2 * DW(N)
;nd0: 1
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: Echo-AntiEcho


;use gradient ratio:	gp 0 : gp 1 : gp 2
;			-16 :  11 :    7


;for z-only gradients:
;gpz0: -16%
;gpz1: 11%
;gpz2:  7%
;gpz3: -23%

;use gradient files:   
;gpnam0: SMSQ10.100
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100




;Processing

;PHC0(F1): 90
;PHC1(F1): -180
;FCOR(F1): 1




