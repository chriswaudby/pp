;1H CEST (starting on zz) for isolated methyl groups
;13C evolution, non-CT
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

;list of CEST saturation frequencies
;give values in Hz relative to sfo1
define list<frequency> H1sat = <$FQ1LIST>


"p2=p1*2"
"d2=p2"
"p4=p3*2"
"d4=1s/(cnst2*4)"
"d11=30m"

"DELTA1=d4-p16-d16-4u-0.5*larger(p2,p4)-0.6366*p1"
"DELTA2=d4-p16-d16-8u-0.5*larger(p2,p4)"

"in0=inf2/2"
"d0=in0/2-p1-p3*0.6366"

"acqt0=0"
baseopt_echo

aqseq 312

1 ze 
  d11 pl1:f1 pl12:f2
2 d11 do:f2
  4u BLKGRAD

# ifdef OFFRES_PRESAT
    30u fq=cnst21(bf hz):f1
# endif /*OFFRES_PRESAT*/
  20u pl9:f1
  d1 cw:f1 ph29
  4u do:f1
  20u pl1:f1 pl2:f2
  30u fq=0:f1
  50u UNBLKGRAD

  (p3 ph1):f2    ; crush eq'm magnetisation
  20u
  p16:gp1
  d16

  ; CEST period
  4u pl8:f1
  4u H1sat:f1
  4u LOCKH_ON
  d18 cw:f1 ph1
  4u do:f1
  4u LOCKH_OFF
  4u pl1:f1 pl2:f2
  4u fq=0:f1
  4u UNBLKGRAD

  ; t1 evolution
  (p3 ph11):f2
  d0
  (p2 ph1) ; central 1H refocusing
  d0
  (p3 ph12):f2

  ; zz crusher
  20u
  p16:gp2
  d16

  ; back transfer
  (p1 ph1) 
  4u
  p16:gp3
  d16
  DELTA1
  (center (p2 ph2) (p4 ph1):f2 )
  DELTA2
  4u pl12:f2
  p16:gp3
  d16 
  4u BLKGRAD

  go=2 ph31 cpd2:f2 
  d11 do:f2 mc #0 to 2
    F1QF(H1sat.inc)
    F2PH(ip11, id0)

exit 
  

ph1=0
ph2=1
ph11=0 2
ph12=0 0 2 2
ph29=0
ph31=0 2 2 0





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



