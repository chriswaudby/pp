;Mar 2017: added CPD and baseopt
;
;Sep 2014: added option for alternative gradient ramp file
;
;May 2013: separate power levels for water flipback and flipdown pulses
;Adjusted delays for zero first-order phase correction
;
;With possibility for multiple acquistion blocks when NS > phase cycle
;  => set TD0 > 1 (total scans = TD0*NS)
;
;From MH_XSte
;Modified to use convention that d20 is equal to big delta
;Reduced time between gradient pulses in bipolar pairs (tau):
;    tau = d16 + p22
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

# ifdef ALT_RAMP
define list<gradient> diff=<Difframp2>
# else
define list<gradient> diff=<Difframp>
# endif /*ALT_RAMP*/

"p2=p1*2"
"p22=p21*2"

"d11=30m"
"d12=20u"
"d13=4u"
"d15=50u"
;"d4=1s/(cnst4*4)-p30-d16-larger(p21,p1)"
;"d5=1s/(cnst4*4)-p19-d16-larger(p21,p1)"
;"d6=1s/(cnst4*4)-p19-d16-larger(p21,p1)-p11-d12"
"d4=2.77m-p30-d16-larger(p21,p1)"
"d5=2.77m-p19-d16-larger(p21,p1)"
"d6=2.77m-p19-d16-larger(p21,p1)-p11-d12"

"DELTA1=d20-8*d16-6*p19-4*p21-3*larger(p22,p2)-3*d5-2*p11-2*d15-2*p30-2*p1-2*d12-2*d4-d6-d13"

"TAU=p1*0.63662+de"
"acqt0=de"

1 ze
  d11
  d12 BLKGRAD
2 d11 do:f3
3 d11
4 d11

# ifdef CRUSHER

  50u UNBLKGRAD
  p19:gp0
  d16
  10u pl1:f1
  (p1 ph4):f1				; +x
  4u pl0:f1
  (p11:sp11 ph9:r):f1			; flipback(-x): -y -> +z
  4u
  p19:gp0*0.71
  d16
  4u BLKGRAD

# endif /*CRUSHER*/

# ifdef PRESAT
  d12 pl9:f1
  d1 cw:f1 ph29
  d13 do:f1
# else
  d1
# endif /*PRESAT*/

  d12 pl3:f3
  (p21 ph4):f3
  d15 UNBLKGRAD
  p19:gp4                     		;Eqm Nz spoiler
  d16 pl0:f1
  (p11:sp1 ph2:r):f1		; flipdown(-x): +z -> +y
  d12 pl1:f1
  (p1 ph4):f1
  d4
  p30:gp6*diff                		;gradient encoding
  d16
  (center (p2 ph4):f1 (p22 ph4):f3)
  p30:gp6*-1*diff                     ;gradient encoding
  d16
  d4
  (p1 ph1):f1
  d12 pl0:f1
  (p11:sp1 ph3:r):f1		; flipdown(-y): -x -> -z
  d15
  p19:gp2                     		;2HzNz spoiler
  d16
  (p21 ph4):f3
  d5
  p19:gp8                     		;spoiler (echo)
  d16 pl1:f1
  (center (p2 ph2):f1 (p22 ph4):f3)
  d5
  p19:gp8                     		;spoiler (echo)
  d16
  (p21 ph5):f3
  d15
  p19:gp3                     		;Nz spoiler
  d16 BLKGRAD
  DELTA1
  (p21 ph6):f3
  d5 UNBLKGRAD
  p19:gp9                     		;spoiler (echo)
  d16
  (center (p2 ph4):f1 (p22 ph4):f3)
  d6
  p19:gp9                     		;spoiler (echo)
  d16 pl0:f1
  (p11:sp11 ph4:r):f1		; flipback(+x): -z -> +y
  d12 pl1:f1
  (p21 ph2):f3
  d13
  (p1 ph2):f1
  d4
  p30:gp6*diff                     	;gradient decoding
  d16
  (center (p2 ph2):f1 (p22 ph4):f3)
  p30:gp6*-1*diff                     ;gradient decoding
  d16
  d4

; Watergate detection

  TAU
  10u pl18:f1
  p16:gp1
  d16
  p27*0.231 ph7
  d19*2
  p27*0.692 ph7
  d19*2
  p27*1.462 ph7
  d19*2
  p27*1.462 ph8
  d19*2
  p27*0.692 ph8
  d19*2
  p0*0.231 ph8
  6u
  p16:gp1
  d16 pl16:f3

  4u BLKGRAD
  go=2 ph31 cpd3:f3
  4u do:f3
  d11 wr #0 if #0 zd igrad diff
  lo to 3 times td1
  d11 rf #0
  lo to 4 times td0
exit


ph1= 1
ph2= 2
ph3= 3
ph4= 0
ph5= 1 3
ph6= 1 1 3 3
ph7= 3 3 3 3 0 0 0 0 1 1 1 1 2 2 2 2
ph8= 1 1 1 1 2 2 2 2 3 3 3 3 0 0 0 0
ph9= 2
ph29=0
ph31=0 2 2 0 2 0 0 2 0 2 2 0 2 0 0 2


;pl1 : f1 channel - power level for pulse (default)
;pl3 : f3 channel - power level for pulse (default)
;pl9 : f1 channel - power level for presaturation
;pl16: f3 channel - power level for CPD/BB decoupling
;pl18: f1 channel - high power pulse, fine-tuning for watergate
;p1 : f1 channel -  90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;p21 : f3 channel -  90 degree high power pulse
;p22 : f3 channel - 180 degree high power pulse
;p11 : f1 channel - water flipdown/back [1-2 ms]
;p27: f1 channel -  90 degree high power pulse, fine-tuning for watergate
;p0 : f1 channel -  90 degree high power pulse, fine-tuning for watergate
;p30 : encode/decode gradient [small delta = 2*p30]
;p16 : watergate gradient
;p19 : spoiler gradient [500-1000 us]
;sp1 : f1 channel - power level for water flipdown
;sp11 : f1 channel - power level for water flipback
;spnam1 : sinc1.1000
;spnam11 : sinc1.1000
;d1 : relaxation delay; 1-5 * T1
;d16 : gradient recovery delay [200usec]
;d19 : delay for binomial water suppression
;     d19 = (1/(2*d)), d = distance of next null (in Hz)
;d20 : diffusion delay (big delta)
;cpd3: decoupling according to sequence defined by cpdprg3
;pcpd3: f3 channel - 90 degree pulse for decoupling sequence
;cnst4: = J(NH)
;NS: 4 * n
;DS: 8
;td0: dimension of accumulation loop (no. of acqusition blocks)
;td1: number of experiments

;gpz0: d1 crusher [73.73 %]
;gpz1: Watergate [-53 %]
;gpz2: 2HzNz crush [17 %]
;gpz3: Nz crush [13 %]
;gpz4: Eqm Nz crush [11 %]
;gpz6: Diffusion [100 %]
;gpz8: 180 pair [9 %]
;gpz9: 180 pair [15 %]

;gpnam0: SMSQ10.100
;gpnam1: SMSQ10.100
;gpnam2: SMSQ10.100
;gpnam3: SMSQ10.100
;gpnam4: SMSQ10.100
;gpnam6: SMSQ10.100
;gpnam8: SMSQ10.100
;gpnam9: SMSQ10.100
