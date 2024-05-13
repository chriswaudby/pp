;mlevgpph_watercontrol.cw
;avance-version (10/01/28)
;homonuclear Hartman-Hahn transfer using MLEV17 sequence for mixing
;using two power levels for excitation and spinlock
;phase sensitive
;water suppression using watergate W5 pulse sequence with gradients
;using double echo
;
;A. Bax & D.G. Davis, J. Magn. Reson. 65, 355-360 (1985)
;M. Liu, X. Mao, C. He, H. Huang, J.K. Nicholson & J.C. Lindon,
;   J. Magn. Reson. 132, 125 - 129 (1998)
;
;$CLASS=HighRes
;$DIM=2D
;$TYPE=
;$SUBTYPE=
;$COMMENT=


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p5=p6*.667"
"p7=p6*2"
"d12=20u"


"in0=inf1"

"d0=in0/2-p1*2/3.1416-4u"


"SCALEF=p7*2/p5"
"FACTOR1=((d9-p17*2)/(p6*64+p5))/SCALEF"
"l1=FACTOR1*SCALEF"


1 ze 
2 d1 
3 d12 pl1:f1
  p1 ph1
  d0
  4u pl10:f1
  (p17 ph26)
						;begin MLEV17
4 (p6 ph22 p7 ph23 p6 ph22)
  (p6 ph24 p7 ph25 p6 ph24)
  (p6 ph24 p7 ph25 p6 ph24)
  (p6 ph22 p7 ph23 p6 ph22)
  (p6 ph24 p7 ph25 p6 ph24)
  (p6 ph24 p7 ph25 p6 ph24)
  (p6 ph22 p7 ph23 p6 ph22)
  (p6 ph22 p7 ph23 p6 ph22)
  (p6 ph24 p7 ph25 p6 ph24)
  (p6 ph22 p7 ph23 p6 ph22)
  (p6 ph22 p7 ph23 p6 ph22)
  (p6 ph24 p7 ph25 p6 ph24)
  (p6 ph22 p7 ph23 p6 ph22)
  (p6 ph22 p7 ph23 p6 ph22)
  (p6 ph24 p7 ph25 p6 ph24)
  (p6 ph24 p7 ph25 p6 ph24)
  (p5 ph23)
  lo to 4 times l1 
						;end MLEV17
  (p17 ph26) 

 50u pl18:f1 UNBLKGRAD
  p16:gp1
  d16
  p27*0.087 ph6
  d19*2
  p27*0.206 ph6
  d19*2
  p27*0.413 ph6
  d19*2
  p27*0.778 ph6
  d19*2
  p27*1.491 ph6
  d19*2
  p27*1.491 ph7
  d19*2
  p27*0.778 ph7
  d19*2
  p27*0.413 ph7
  d19*2
  p27*0.206 ph7
  d19*2
  p27*0.087 ph7
   p16:gp2
  50u BLKGRAD
  d16  
p1 ph10
  50u UNBLKGRAD 
  p17:gp3
  50u BLKGRAD
  d16
  d25
p1 ph11
   50u UNBLKGRAD 
   p16:gp1
  d16
  p27*0.087 ph8
  d19*2
  p27*0.206 ph8
  d19*2
  p27*0.413 ph8
  d19*2
  p27*0.778 ph8
  d19*2
  p27*1.491 ph8
  d19*2
  p27*1.491 ph9
  d19*2
  p27*0.778 ph9
  d19*2
  p27*0.413 ph9
  d19*2
  p27*0.206 ph9
  d19*2
  p27*0.087 ph9
   p16:gp2
  50u BLKGRAD
  d16 
  go=2 ph31
  d1 mc #0 to 2 F1PH(ip1,id0)
exit

  
ph1=0 2
ph6=0 0 2 2
ph7=2 2 0 0
ph8=0 0 2 2
ph9=2 2 0 0
ph10=1
ph11=3
ph22=3
ph23=0
ph24=1
ph25=2
ph26=0
ph31=0 2 


;pl1 : f1 channel - power level for pulse (default)
;pl10: f1 channel - power level for TOCSY-spinlock
;pl18: f1 channel - power level for W5-pulse (watergate)
;p1 : f1 channel -  90 degree high power pulse
;p5 : f1 channel -  60 degree low power pulse
;p6 : f1 channel -  90 degree low power pulse
;p7 : f1 channel - 180 degree low power pulse
;p16: homospoil/gradient pulse
;p17: f1 channel -  trim pulse                       [2.5 msec]
;p27: f1 channel -  90 degree pulse at pl18
;d0 : incremented delay (2D)
;d1 : relaxation delay; 1-5 * T1
;d9 : TOCSY mixing time
;d12: delay for power switching                      [20 usec]
;d16: delay for homospoil/gradient recovery
;d19: delay for binomial water suppression
;     d19 = (1/(2*d)), d = distance of next null (in Hz)
;l1: loop for MLEV cycle: (((p6*64) + p5) * l1) + (p17*2) = mixing time
;inf1: 1/SW = 2 * DW
;in0: 1/(1 * SW) = 2 * DW
;nd0: 1
;NS: 2 * n
;DS: 16
;td1: number of experiments
;FnMODE: States-TPPI, TPPI, States or QSEQ

;use gradient ratio:	gp 1 : gp 2
;			  34 :   22

;for z-only gradients:
;gpz1: 34%
;gpz2: 22%

;use gradient files:   
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100


;Processing

;PHC0(F1): 180
;PHC1(F1): -180
;FCOR(F1): 1



;$Id: mlevgpphw5,v 1.9.2.1 2010/02/02 15:31:39 ber Exp $
