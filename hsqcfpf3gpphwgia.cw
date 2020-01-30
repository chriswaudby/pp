;Modified to use half-dwell initial sampling delay by default (May 2013)
;
;With in-phase/anti-phase element, JK, July 2012
;Crushing equm 15N magnetisation
;
;Without refocusing for 0,0 phase correction
;Option for either (90,-180) or (180,-360) phase correction
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
"d26=1s/(cnst4*4)"


"in0=inf1/4"

"DELTA1=d26-p16-d16"
"DELTA2=d26-p19-d16-p11-12u"
"DELTA5=d26-p19-d16-p11-12u-de+0.63662*p1"

"DELTA3=d26/2-p19-d16-10u"
"DELTA4=d26-p19-d16-10u"

# ifdef SINGLEDWELL
    "d0=in0-0.5*(10u+0.63622*p21)"
# else
    "d0=in0/2-0.5*(10u+0.63622*p21)"
# endif /*SINGLEDWELL*/

"l0=1"


1 ze 
  d11 pl16:f3
2 d11 do:f3

# ifdef LABEL_CN
  if "4*d0-24u > p8"
  {
   d12  
   "d60=d0-0.25*p8"
  }
  d12 pl0:f2
# endif /*LABEL_CN*/

  d1

3 d12 pl1:f1 pl3:f3
  50u UNBLKGRAD
  (p21 ph1):f3
  p16:gp1
  d16
  (p21 ph2):f3
  p16:gp1*0.7
  d16*2

  (p1 ph1)
  p16:gp2
  d16
  DELTA1
  (center (p2 ph2) (p22 ph6):f3 )
  DELTA1
  p16:gp2
  d16
  (p1 ph2) 

  4u pl0:f1
  (p11:sp11 ph8:r):f1	; flipback(-x), -y -> +z
  4u
  p16:gp3
  d16 pl1:f1

; in-phase/anti-phase element

  if "l0 %2 == 1"      ; in-phase (refocusing)
    {
    (p21 ph10):f3
    d26*0.5
    (p2 ph5):f1
    10u
    DELTA3
    p19:gp4
    d16
    (p22 ph1):f3
    10u
    p19:gp4
    d16
    DELTA3
    (p2 ph5):f1
    d26*0.5
    (p21 ph10):f3
    }
  else                  ; anti-phase (evolving)
    {
    (p21 ph10):f3
    DELTA4
    10u
    p19:gp4
    d16
    (p2 ph5):f1
    2u
    (p22 ph1):f3
    10u
    p19:gp4
    d16
    DELTA4
    (p2 ph5):f1
    2u
    (p21 ph11):f3
    }

; end of IPAP element

  10u
  p16:gp5
  d16

# ifdef LABEL_CN

  if "4*d0-24u < p8"
  { 
   (p21 ph3):f3
   2u
   d0 gron0
   d0 gron0*-1
   8u groff
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
   (p8:sp13 ph1):f2
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
  2u
  d0 gron0
  d0 gron0*-1
  8u groff
  (p21 ph4):f3

# endif /*LABEL_CN*/

  4u
  p16:gp6
  d16 pl0:f1
  (p11:sp1 ph7:r):f1	; flipdown(-x), z -> y
  4u
  4u pl1:f1

  (p1 ph1) 
  4u
  p19:gp7
  d16
  DELTA2 pl0:f1
  (p11:sp1 ph7:r):f1	; flipdown(-x), z -> y
  4u
  4u pl1:f1
  (center (p2 ph1) (p22 ph1):f3 )
  4u pl0:f1
  (p11:sp11 ph9:r):f1	; flipback(-x), -y -> z
  4u
  p19:gp7
  d16
  DELTA5 pl16:f3
  4u BLKGRAD

  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 
	F1I(iu0, 2)
	F1PH(ip3 & ip6, id0)
exit 
  

ph1=0
ph2=1
ph3=0 2
ph4=0 0 2 2
ph5=0 0 0 0 2 2 2 2
ph6=0
ph7=2
ph8=2
ph9=2
ph10=0
ph11=1
ph31=0 2 2 0


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
;d16: delay for homospoil/gradient recovery
;d26 : 1/(4J)YH
;cnst4: = J(YH)
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
;gpz1: 47%
;gpz2: 13%
;gpz3: 31%
;gpz4: 17%
;gpz5: 41%
;gpz6: 47%
;gpz7: 53%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam5: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam7: SMSQ10.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: hsqcfpf3gpphwg,v 1.6.2.1 2007/07/04 13:41:19 ber Exp $
