;best_trosy_bshd.cw
;
;based on b_trosyf3gpph.2 / homodec1D_dqd_mlev.jfy
;avance-version (15/03/12)
;best-TROSY
;2D H-1/X correlation via TROSY
;phase sensitive using Echo/Antiecho method
;using f3 - channel
;with additional 180degree pulse on N-15
;(use parameterset B_TROSYF3GPPH)
;
;A. Favier & B. Brutscher, J. Biomol. NMR 49, 9-15 (2011)
;E. Lescop, T. Kern & B. Brutscher, J. Magn. Reson. 203, 190-198 (2010)
;(P. Schanda, H. v. Melckebeke & B. Brutscher,
;   J. Am. Chem. Soc. 128, 9042-9043 (2006))
;(E. Lescop, P. Schanda & B. Brutscher,
;   J. Magn. Reson.  187 163-169 (2007))
;(M. Czisch & R. Boelens, J. Magn. Reson. 134, 158-160 (1998))
;(K. Pervushin, G. Wider & K. Wuethrich, J. Biomol. NMR 12,
;   345-348 (1998))
;(A. Meissner, T. Schulte-Herbrueggen, J. Briand & O.W. Sorensen, Mol. Phys. 96,
;   1137-1142 (1998))
;(J. Weigelt, J. Am. Chem. Soc. 120, 10778-10779 (1998))
;(M. Rance, J.P. Loria & A.G. Palmer III, J. Magn. Reson. 136, 91-101 (1999))
;(G. Zhu, X.M. Kong & K.H. Sze, J. Biomol. NMR 13, 77-81 (1999))
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
#include <De.incl>


; bshd bits
; NB acquisition time must be a multiple of 48 ms
;#define BASH
define loopcounter tdCount

"d27=12m"   ;half length between the two decoupling blocks
"l27=d27/(dw*2)+0.5"

"d28=dw*2*l27"
"d29=d28*2"

"tdCount = aq/(d28*4)+1"

"d3=de+12u"

"l24=0"
"cnst22=0"
"cnst24=0"

dwellmode explicit




; trosy bits
"p22=p21*2"
"d11=30m"
"d12=20u"
"d26=2.7m"

"p29=250u"


"d0=3u"
"in0=inf1/2"


"DELTA1=d26-p19-d16-larger(p22,p42)/2"
"DELTA2=d26-p29-d16-larger(p22,p42)/2-p43*cnst43"
"DELTA3=d26-p16-d16-larger(p22,p42)/2"
"DELTA4=de"

#   ifdef LABEL_CN
"DELTA=d0*2+p8+p21*4/PI"
#   else
"DELTA=d0*2+p21*4/PI"
#   endif /*LABEL_CN*/


"l0=1"


"spoff13=bf2*(cnst26/1000000)-o2"
"spoff26=bf1*(cnst54/1000000)-o1"
"spoff28=bf1*(cnst54/1000000)-o1"
"spoff29=bf1*(cnst54/1000000)-o1"



1 d11 ze
2 d11
3 d12

  (p22 ph1):f3
  d1 pl9:f1
  d18 cw:f1 ph9
  10u do:f1 pl1:f1
  10u UNBLKGRAD
  p16:gp0
  d16

  (p43:sp28 ph1)
  p19:gp1
  d16
  DELTA1
  (center (p42:sp26 ph1) (p22 ph1):f3 )
  DELTA1
  p19:gp1
  d16
  (p43:sp29 ph2)

  p16:gp2
  d16

  if "l0 %2 == 1"
     {
     (p21 ph5):f3
     }
  else
     {
     (p21 ph6):f3
     }

  d0

#   ifdef LABEL_CN
  (p8:sp13 ph1):f2
#   else
#   endif /*LABEL_CN*/

  d0
  (p22 ph1):f3
  DELTA

  (p43:sp29 ph7)
  p29:gp4
  d16
  DELTA2
  (center (p42:sp26 ph1) (p22 ph1):f3 )
  DELTA2
  p29:gp4
  d16
  (p43:sp28 ph1)

  DELTA4
  (p21 ph7):f3
  4u
  p16:gp3
  d16
  DELTA3
  (center (p42:sp26 ph1) (p22 ph1):f3 )
  DELTA3
  p16:gp3
  d16
  4u BLKGRAD
  (p21 ph1):f3

  ; acquisition with BASHD
  ACQ_START(ph30,ph31) ;total delay here=de

   0.05u DWL_CLK_ON
   0.1u REC_UNBLK
77  d28

#ifdef BASH

   if "l24%8==0 || l24%8==5"
      {
       "cnst22=0"
       "cnst24=0"
      }
   if "l24%8==1 || l24%8==4"
      {
       "cnst22=180"
       "cnst24=180"
      }
   if "l24%8==3 || l24%8==6"
      {
       "cnst22=0"
       "cnst24=180"
      }
   if "l24%8==2 || l24%8==7"
      {
       "cnst22=180"
       "cnst24=0"
      }

   0.1u REC_BLK
   0.05u DWL_CLK_OFF

   5u
   p31:gp6
   20u ip12+cnst22
   23u ip14+cnst22 ;may need a longer delay, depending on your grad recovery
   3u pl0:f1
  (p10:sp0 ph12):f1
   6u
   p31:gp6
   45u            ;may need a longer delay, depending on your grad recovery
   p32:gp7
   40u            ;may need a longer delay, depending on your grad recovery
   60u pl1:f1
  (p1*2 ph14):f1
   5u
   p32:gp7
   95u            ;may need a longer delay, depending on your grad recovery
   0.05u DWL_CLK_ON
   0.1u REC_UNBLK
#endif

   d29

#ifdef BASH
   0.1u REC_BLK
   0.05u DWL_CLK_OFF

   5u
   p31:gp8
   20u ip12+cnst24
   20u ip14+cnst24 ;may need a longer delay, depending on your grad recovery
   5u pl0:f1
  (p10:sp0 ph12):f1
   5u
   p31:gp8
   45u            ;may need a longer delay, depending on your grad recovery
   p32:gp9
   40u            ;may need a longer delay, depending on your grad recovery
   60u pl1:f1
  (p1*2 ph14):f1
   5u
   p32:gp9
   95u iu24       ;may need a longer delay, depending on your grad recovery
   0.05u DWL_CLK_ON
   0.1u REC_UNBLK
#endif

   d28

lo to 77 times tdCount

   0.1u REC_BLK
   0.05u DWL_CLK_OFF

   6u ru24 

  rcyc = 2
;  50m wr #0 
  d11 mc #0 to 2 
     F1EA(calph(ph7, +180) & calclc(l0, 1), caldel(d0, +in0) & calph(ph5, +180) & calph(ph6, +180) & calph(ph31, +180))
exit


ph1=0
ph2=1 
ph3=2
ph4=3
ph5=1 3 0 2
ph6=1 3 2 0
ph7=3
ph9=0
ph12=2
ph14=0
ph30=0
ph31=1 3 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;sp13: f2 channel - shaped pulse 180 degree (Ca and C=O, adiabatic)
;sp26: f1 channel - shaped pulse 180 degree (Reburp.1000)
;sp28: f1 channel - shaped pulse  90 degree (Eburp2.1000)
;sp29: f1 channel - shaped pulse  90 degree (Eburp2tr.1000)
;                   for time reversed pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse                         [1 msec]
;p19: gradient pulse 2                                 [500 usec]
;p31: gradient 1 for bashd  [210 usec]
;p32: gradient 2 for bashd  [210 usec]
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p29: gradient pulse 3                                 [250 usec]
;p42: f1 channel - 180 degree shaped pulse for refocussing
;                      Reburp.1000               (1.4ms at 600.13 MHz)
;p43: f1 channel -  90 degree shaped pulse for excitation
;                      Eburp2.1000/Eburp2tr.1000 (1.7ms at 600.13 MHz)
;d0 : incremented delay (F1)                           [3 usec]
;d1 : relaxation delay (excluding presat time)
;d18 : presaturation delay
;d11: delay for disk I/O                               [30 msec]
;d12: delay for power switching                        [20 usec]
;d16: delay for homospoil/gradient recovery
;d26: 1/(4J(NH)                                        [2.7 msec]
;cnst26: Call chemical shift (offset, in ppm)          [101 ppm]
;cnst43: compensation of chemical shift evolution during p43
;           Eburp2.1000: 0.69
;cnst52: scaling factor for p42 to compensate for transition region
;           Reburp.1000: 1.426
;cnst53: scaling factor for p43 to compensate for transition region
;           Eburp2.1000: 1.000
;cnst54: H(N) chemical shift (offset, in ppm)
;cnst55: H(N) bandwidth (in ppm)
;inf1: 1/SW(N) = 2 * DW(N)
;in0: 1/(2 * SW(N)) = DW(N)
;nd0: 2
;ns: 4 * n
;ds: >= 16
;td1: number of experiments in F1
;FnMODE: echo-antiecho in F1


;for z-only gradients:
;gpz0: 40%
;gpz1: 30%
;gpz2: 60%
;gpz3: 44%
;gpz4: 18%
;gpz6: 14.7%
;gpz7: 21%
;gpz8: 25.9%
;gpz9: 35%

;use gradient files:   
;gpnam0: SMSQ10.100
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.32
;gpnam6: SINE.20
;gpnam7: SINE.20
;gpnam8: SINE.20
;gpnam9: SINE.20


                                          ;preprocessor-flags-start
;CALC_SP: for calculation of all bandselective Proton pulses based on cnst54 and cnst55
;             option -DCALC_SP (eda: ZGOPTNS)
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with 
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end
										  


;Processing

;PHC0(F1): -22.5
										  


;$Id: b_trosyf3gpph.2,v 1.2.2.3 2015/03/12 17:07:08 ber Exp $

