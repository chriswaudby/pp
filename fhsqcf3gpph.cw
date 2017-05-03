;fhsqcf3gpph
; Crushing equm 15N magnetisation
; Delays adjusted for (90,180) phase correction
; option for 13C decoupling during acquisition (-DDUAL_DECOUPLING)
;
;avance-version (07/04/04)
;2D H-1/X correlation via double inept transfer
;phase sensitive
;with decoupling during acquisition
;
;S. Mori, C. Abeygunawardana, M. O'Neil-Johnson & P.C.M. van Zijl,
;   J. Magn. Reson. B 108, 94-98 (1995)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
"p22=p21*2"
"d11=30m"
"d12=20u"
"d13=4u"
"d26=1s/(cnst4*4)"


"in0=inf1/2"

"d0=in0-0.5*larger(p2,p14)-0.6366*p21"


"DELTA=d19-p22/2"
"DELTA1=d26-p16-d16-p27*3-d19*5-p1*2/PI"
"DELTA2=d26-p16-d16-p27*2-p0-d19*5-de-8u"
"DELTA4=p21*2/3.1416"
# ifdef LABEL_CN
"DELTA3=d0+larger(p2,p14)/2"
# else
"DELTA3=d0+p2/2"
# endif /*LABEL_CN*/


"TAU=d26-p16-4u"


"acqt0=0"
baseopt_echo


1 ze
  d11 pl16:f3 pl12:f2
2 d11 do:f3 do:f2

  d1

3 d12 pl1:f1 pl3:f3
  50u UNBLKGRAD

  ; crush equilibrium 15N
  p16:gp4*0.7
  d16
  (p21 ph1):f3
  p16:gp4
  d16

  ; begin HSQC
  (p1 ph1)
  4u
  p16:gp1
  TAU pl3:f3
  (center (p2 ph1) (p22 ph6):f3 )
  4u
  p16:gp1
  TAU
  (p1 ph2)

  ; zz filter
  4u
  p16:gp2
  d16

  ; t1 evolution

  (p21 ph3):f3
;  DELTA3
;  (p22 ph3):f3
;  DELTA4
  d0

#   ifdef LABEL_CN
  (center (p2 ph5) (p14:sp3 ph1):f2 )
#   else
  (p2 ph5)
#   endif /*LABEL_CN*/

  d0
;  DELTA4
;  (p22 ph4):f3
;  DELTA3
  (p21 ph4):f3

  ; zz filter
  4u
  p16:gp2
  d16

  ; back transfer
  (p1 ph7)
  DELTA1
  p16:gp3
  d16 pl18:f1
  p27*0.231 ph8
  d19*2
  p27*0.692 ph8
  d19*2
  p27*1.462 ph8
  DELTA
  (p22 ph1):f3
  DELTA
  p27*1.462 ph9
  d19*2
  p27*0.692 ph9
  d19*2
  p0*0.231 ph9
  4u
  p16:gp3
  d16
  4u BLKGRAD

  ; acquitision with dual decoupling
# ifdef DUAL_DECOUPLING
  DELTA2 pl12:f2 pl16:f3
  go=2 ph31 cpd2:f2 cpd3:f3
  d11 do:f2 do:f3 mc #0 to 2 F1PH(ip3 & ip6, id0)
# else
  DELTA2 pl16:f3
  go=2 ph31 cpd3:f3
  d11 do:f3 mc #0 to 2 F1PH(ip3 & ip6, id0)
# endif /*LABEL_CN*/

exit


ph1=0
ph2=1
ph3=0 2
ph4=0 0 0 0 2 2 2 2
ph5=0 0 2 2
ph6=0
ph7=2
ph8=1
ph9=3
ph31=0 2 0 2 2 0 2 0


;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl16: f3 channel - power level for CPD/BB decoupling
;pl18: f1 channel - power level for 3-9-19-pulse (watergate)
;sp3: f2 channel - shaped pulse 180 degree (adiabatic)
;p0 : f1 channel -  90 degree pulse at pl18
;                      use for fine adjustment
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p14: f2 channel - 180 degree shaped pulse for inversion (adiabatic)
;p16: homospoil/gradient pulse
;p21: f3 channel -  90 degree high power pulse
;p22: f3 channel - 180 degree high power pulse
;p27: f1 channel -  90 degree pulse at pl18
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O                             [30 msec]
;d12: delay for power switching                      [20 usec]
;d13: short delay                                    [4 usec]
;d16: delay for homospoil/gradient recovery
;d19: delay for binomial water suppression
;     d19 = (1/(2*d)), d = distance of next null (in Hz)
;d26 : 1/(4J(YH))
;cnst4: = J(YH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 8 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI (or TPPI)
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence


;for z-only gradients:
;gpz1: 19%
;gpz2: 73%
;gpz3: 31%
;gpz4: 53%

;use gradient files:
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SINE.100


                                          ;preprocessor-flags-start
;LABEL_CN: for C-13 and N-15 labeled samples start experiment with
;             option -DLABEL_CN (eda: ZGOPTNS)
                                          ;preprocessor-flags-end


;$Id: fhsqcf3gpph,v 1.6 2007/04/11 13:34:29 ber Exp $
