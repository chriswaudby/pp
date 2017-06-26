;With 13C decoupling: check AQ < 140 ms
;
;With possibility for multiple acquisition blocks when NS > phase cycle
;  => set TD0 > 1 (total scans = TD0*NS)
;
;From stebpgp1s19xc.jk
;Using presaturation for water suppression
;Adiabatic inversion pulses on 13C during INEPT periods where proton magnetization transverse.
;Hard 180s on 13C during INEPT periods where carbon magnetization transverse.
;Moved water flipback into ZZ period during back transfer
;
;From MH_XSte
;Modified to use convention that d20 is equal to big delta
;Reduced time between gradient pulses in bipolar pairs (tau):
;    tau = d16 + p8
;
;H-1/X correlation via double refocused inept transfer
;ste during the transfer steps and storage of the magnetization on the X-nucleus during the diffusion delay
;watergate after the decoding gradients for use with z-only gradient probes
;1D version
;written by Fabien Ferrage, last modification November 22nd 2004
;
;Ferrage et al., JACS (2004) 126:5654

prosol relations=<triple_d>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<gradient> diff=<Difframp>

"p2=p1*2"
"p4=p3*2"
"d4=1s/(cnst4*4)-p30-d16-0.5*larger(p8,p2)"
"d5=1s/(cnst4*cnst12)-p19-d16-0.5*larger(p4,p2)"
"d11=30m"
"d13=4u"
"d12=20u"
"d15=50u"

"DELTA1=d20-8*d16-6*p19-4*p3-2*larger(p4,p2)-4*d5-2*d12-2*d13-2*d15-2*p30-2*p1-2*d4-p8"

1 ze
  d11
  d12 BLKGRAD
2 d11 do:f2
3 d11
4 d12
  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1 pl2:f2
  (p3 ph4):f2
  d15 UNBLKGRAD
  p19:gp4                     		;Eqm Cz spoiler
  d16
  (p1 ph4):f1
  d4
  p30:gp6*diff                		;gradient encoding
  d16 pl0:f2
  (center (p2 ph4):f1 (p8:sp13 ph4):f2)
  p30:gp6*-1*diff                     ;gradient encoding
  d16 pl2:f2
  d4
  (p1 ph1):f1
  d15
  p19:gp2                     		;2HzCz spoiler
  d16
  (p3 ph4):f2
  d5
  p19:gp8                     		;spoiler (echo)
  d16 
  (center (p2 ph2):f1 (p4 ph4):f2)
  d5
  p19:gp8                     		;spoiler (echo)
  d16
  (p3 ph5):f2
  d15
  p19:gp3                     		;Cz spoiler
  d16 BLKGRAD
  d12 pl9:f1
  DELTA1 cw:f1 ph29
  d13 do:f1
  d12 pl1:f1
  (p3 ph6):f2
  d5 UNBLKGRAD
  p19:gp9                     		;spoiler (echo)
  d16
  (center (p2 ph4):f1 (p4 ph4):f2)
  d5 
  p19:gp9                     		;spoiler (echo)
  d16 
  (p3 ph2):f2
  d13
  (p1 ph2):f1 
  d4
  p30:gp6*diff                     	;gradient decoding
  d16 pl0:f2
  (center (p2 ph2):f1 (p8:sp13 ph4):f2)
  p30:gp6*-1*diff                     ;gradient decoding
  d16 pl12:f2
  d4 BLKGRAD
  go=2 ph31 cpd2:f2
  d11 do:f2 wr #0 if #0 zd igrad diff
  lo to 3 times td1
  d11 do:f2 rf #0
  lo to 4 times td0
exit


ph1= 1
ph2= 2
ph3= 3
ph4= 0
ph5= 1 3
ph6= 1 1 3 3
ph29=0
ph31=0 2 2 0 


;pl1 : f1 channel - power level for pulse (default)
;pl2 : f2 channel - power level for pulse (default)
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p3 : f2 channel -  90 degree high power pulse
;p4 : f2 channel - 180 degree high power pulse
;p8 : f2 channel - adiabatic inversion pulse
;p30 : encode/decode gradient [small delta = 2*p30]
;p19 : spoiler gradient [500-900 us]
;sp13: f2 channel - power level for adiabatic inversion
;spnam13 : Crp60,0.5,20.1
;d1 : relaxation delay; 1-5 * T1
;d16 : gradient recovery delay [200usec]
;d20 : diffusion delay (big delta)
;cnst4: = J(CH) [125-145 Hz]
;cnst12: multiplicity selection - 4=CH,CH3; 8=CH2,CH3,(CH); 10=CH3,CH2,(CH) 
;NS: 4 * n
;DS: 8
;td0: dimension of accumulation loop (no. of acqusition blocks)
;td1: number of experiments

;gpz2: 2HzCz crush [17 %]
;gpz3: Cz crush [13 %]
;gpz4: Eqm Cz crush [11 %]
;gpz6: Diffusion [100%]
;gpz8: 180 pair [9 %]
;gpz9: 180 pair [15 %]

;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam8: SMSQ10.100
;gpnam9: SMSQ10.100
