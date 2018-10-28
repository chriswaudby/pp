;hsqcph.cw
;avance-version (07/04/04)
;HSQC
;2D H-1/X correlation via double inept transfer
;phase sensitive
;with decoupling during acquisition
;
;G. Bodenhausen & D.J. Ruben, Chem. Phys. Lett. 69, 185 (1980)
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
"d2=p2"
"p4=p3*2"
"d4=1s/(cnst2*4)"
"d11=30m"

"DELTA1=d4-p16-d16-4u-0.5*larger(p2,p4)-0.6366*p1"
"DELTA2=d4-p16-d16-8u-0.5*larger(p2,p4)-de"

"in0=inf1/2"
"d0=in0/2-p1-p3*0.6366"

"acqt0=de"
baseopt_echo

1 ze 
  d11 pl1:f1 pl12:f2
2 d1 do:f2

  4u pl2:f2
  50u UNBLKGRAD
  (p3 ph1):f2    ; crush eq'm magnetisation
  20u
  p16:gp1
  d16

  (p1 ph1)
  4u
  p16:gp2
  d16
  DELTA1
  (center (p2 ph2) (p4 ph9):f2 )
  DELTA1
  4u
  p16:gp2
  d16
  (p1 ph3) 

#ifdef CALIB_F2
;crusher
  50u UNBLKGRAD
  p16:gp0*0.9
  d16

  4u pl20:f2
  (p20 ph1):f2
  p16:gp0
  d16 pl2:f2
#endif

  (p3 ph6):f2
  d0
  (p2 ph8) ; central 1H refocusing
  d0
  (p3 ph7):f2

  (p1 ph4) 
  4u
  p16:gp3
  d16
  DELTA1
  (center (p2 ph2) (p4 ph5):f2 )
  DELTA2
  4u pl12:f2
  p16:gp3
  d16 
  4u BLKGRAD

  go=2 ph31 cpd2:f2 
  d1 do:f2 mc #0 to 2 F1PH(ip6 & ip9, id0)

exit 
  

ph1=0
ph2=0
ph3=1
ph4=1
ph5=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph6=0 2
ph7=0 0 0 0 2 2 2 2
ph8=0 0 2 2
ph9=0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2
ph31=0 2 0 2 2 0 2 0





;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;pl12: f2 channel - power level for CPD/BB decoupling
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p16: homospoil/gradient pulse
;d0 : incremented delay (2D)                         [3 usec]
;d1 : relaxation delay; 1-5 * T1
;d4 : 1/(4J)XH
;d11: delay for disk I/O                             [30 msec]
;d16: delay for homospoil/gradient recovery
;cnst2: = J(XH)
;inf1: 1/SW(X) = 2 * DW(X)
;in0: 1/(2 * SW(X)) = DW(X)
;nd0: 2
;NS: 4 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ
;cpd2: decoupling according to sequence defined by cpdprg2
;pcpd2: f2 channel - 90 degree pulse for decoupling sequence



;$Id: hsqcph,v 1.4 2007/04/11 13:34:30 ber Exp $
