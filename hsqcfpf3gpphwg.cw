;Added option for dual 13C/15N decoupling during acquisition, Aug 2014
;
;Added option for offres presat, Aug 2014
;
;Added option for 4-scan phase cycle (-DFOURSCAN), Sep 2013
;
;Added options (May 2013):
;	- coupled spectrum in F1 (-DCOUPLED)
;	- 1D spectrum with no 15N shift evolution (-DONE_D)
;	- calibration of 15N 90deg pulse (-DCAL_N)
;
;Delays adjusted for zero first-order phase correction
;
;Without refocusing for 0,0 phase correction
;Option for (180,-360) phase correction (default is (90,-180) phase corr.)
;Option for carbon decoupling
;Assumes p8 > p2

;With gradients during t1 to keep water along z
;With separate shape/power level for second water pulse in WATERGATE (flip-back)
;Renamed: hsqcfpf3gpphwg.3 --> hsqcfpf3gpphwg.3.jk (19/1/12)
;
;hsqcfpf3gpphwg
;avance-version (07/06/20)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive
;with decoupling during acquisition
;using f3 - channel
;using flip-back pulse
;water suppression using watergate sequence
;similar to fhsqc 
;(use parameterset )
;
;G. Bodenhausen & D.J. Ruben, Chem. Phys. Lett. 69, 185 (1980)
;M. Piotto, V. Saudek & V. Sklenar, J. Biomol. NMR 2, 661 - 666 (1992)
;V. Sklenar, M. Piotto, R. Leppik & V. Saudek, J. Magn. Reson.,
;   Series A 102, 241 -245 (1993)
;S. Mori, C. Abeygunawardana, M. O'Neil-Johnson & P.C.M. van Zijl,
;   J. Magn. Reson. B 108, 94-98 (1995)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


prosol relations=<triple>

  
#include <Avance.incl>
#include <Delay.incl>
#include <Grad.incl>


"p2=p1*2"
"p22=p21*2"

"d11=30m"
"d12=20u"
"d13=4u"
"d26=1s/(cnst4*4)"

"in0=inf1/4"

"DELTA1=d26-p16-d16-4u"
"DELTA2=d26-p19-d16-p11-d12-4u"
"DELTA3=d26-p19-d16-p11-d12-8u-de+0.63662*p1"

# ifndef ONE_D

# ifdef COUPLED
    "d2=p2"
# endif /*COUPLED*/

# ifdef SINGLEDWELL
    "d0=in0-0.5*(10u+p1+0.63662*p21)"
# else
    "d0=in0/2-0.5*(10u+p1+0.63662*p21)"
# endif /*SINGLEDWELL*/

# endif /*ONE_D*/


1 ze 
  d11 pl16:f3

# ifdef LABEL_CN

2 d11 do:f2 do:f3

  if "4*d0+p2-24u > p8"
  {
   d12  
   "d60=d0+0.5*p1-0.25*p8"
   d12 pl0:f2
  }

# else

2 d11 do:f3

# endif /*LABEL_CN*/

# ifdef OFFRES_PRESAT

  30u fq=cnst21(bf hz):f1
  d12 pl9:f1
  d13 cw:f1 ph1
  d1
  d13 do:f1
  30u fq=0:f1
  
# else

  d1

# endif /*OFFRES_PRESAT*/

3 d12 pl1:f1 pl3:f3
 
  50u UNBLKGRAD
  p16:gp1*0.71
  d16
  (p21 ph1):f3
  p16:gp1
  d16*2

  (p1 ph1)
  4u
  p16:gp3
  d16
  DELTA1
  (center (p2 ph2) (p22 ph6):f3 )
  DELTA1
  4u
  p16:gp3
  d16
  (p1 ph2) 

  4u pl0:f1

# ifdef COUPLED
  (p11:sp11 ph9:r):f1	; flipback(-x), -y -> +z
# else
  (p11:sp1 ph8:r):f1	; flipdown(+x), -y -> -z
# endif /*COUPLED*/

  d12 pl1:f1

# ifdef CAL_N

  4u
  p16:gp4*0.3
  d16
  d12 pl23:f3
  (p25 ph1):f3
  4u
  p16:gp4*0.7
  d16 
  d12 pl3:f3

# else

  4u
  p16:gp4
  d16

# endif /*CAL_N*/

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
# ifdef COUPLED
   (d2):f1					; for F1-coupled spectrum
# else
   (p2 ph5):f1
# endif /*COUPLED*/
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
# ifdef COUPLED
   ( center (d2):f1 (p8:sp13 ph1):f2 )		; for F1-coupled spectrum
# else
   ( center (p2 ph5):f1 (p8:sp13 ph1):f2 )
# endif /*COUPLED*/
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
# ifdef COUPLED
  (d2):f1					; for F1-coupled spectrum
# else
  (p2 ph5):f1
# endif /*COUPLED*/
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
  (p11:sp1 ph7:r):f1	; flipdown(-x), z -> y
  4u
  4u pl1:f1

  (p1 ph1)
  4u
  p19:gp5
  d16 pl0:f1
  DELTA2
  (p11:sp1 ph7:r):f1    ; flipdown(-x), z -> y
  d12 pl1:f1
  (center (p2 ph1) (p22 ph1):f3 )
  d12 pl0:f1
  (p11:sp11 ph9:r):f1   ; flipback(-x), -y -> z
  4u

# if defined(LABEL_CN) && defined(DUALDEC)

  p19:gp5
  d16 pl12:f2 pl16:f3
  DELTA3
  4u BLKGRAD
  go=2 ph31 cpd2:f2 cpd3:f3
  d11 do:f2 do:f3 mc #0 to 2 F1PH(ip3 & ip6, id0)

# else

  p19:gp5
  d16 pl16:f3
  DELTA3
  4u BLKGRAD
  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 F1PH(ip3 & ip6, id0)

# endif /*LABEL_CN, DUALDEC*/


exit 
  

ph1=0
ph2=1
ph3=0 2

# ifdef FOURSCAN
ph4=0 0 2 2		; 4-scan phase cycle
# else
ph4=0 0 0 0 2 2 2 2	; 8-scan phase cycle
# endif /*FOURSCAN*/

ph5=0 0 2 2
ph6=0
ph7=2
ph8=0
ph9=2

# ifdef FOURSCAN
ph31=0 2 2 0		; 4-scan phase cycle
# else
ph31=0 2 0 2 2 0 2 0	; 8-scan phase cycle
# endif /*FOURSCAN*/


;pl0 : 120dB
;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl16: f3 channel - power level for CPD/BB decoupling
;sp1: f1 channel - shaped pulse  90 degree (flip-down)
;sp11: f1 channel - shaped pulse  90 degree (flip-back)
;sp13: f2 channel - shaped pulse 180 degree (adiabatic)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p8 : f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p11: f1 channel -  90 degree shaped pulse
;p16: homospoil/gradient pulse
;p19: second homospoil/gradient pulse
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;d26 : 1/(4J)YH
;cnst4: = J(YH)
;cnst21: offset in Hz for off-resonance presaturation
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz0: 1-2%
;gpz1: 41%
;gpz3: 7%
;gpz4: 17%
;gpz5: 53%
;gpz6: 13%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: hsqcfpf3gpphwg,v 1.6.2.1 2007/07/04 13:41:19 ber Exp $
