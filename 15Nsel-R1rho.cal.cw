;15Nsel-R1rho_HC1.alh
;fred 8-2011
;Selective 1D experiment for measuring N-15 R1rho relaxation times
;
; Hansen, DF, Vallurupalli, P, Kay, LE. J PChem B 112, 5898-5904, (2008) 
; Korzhnev, DM, Orekhov, VY, Kay, LE. JACS 127, 713 (2005)

;phaseOnTcu

prosol relations=<triple>


#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


"p2=p1*2"
;"p4=p3*2"
;"p22=p21*2"

"d0=5u"
"d11=30m"
"d12=20u"
"p31=d31"

"cnst29=0"
"nbl=td1"

#   ifdef ZETA_DELAY
"d2=1s/(4*cnst10)"
#	else
"d2=0s"
#   endif /*ZETA_DELAY*/ 

"p29=d30-p31"
"d29=d30-p31"

"p18=p1*atan(cnst11/(cnst30*bf1-(sfo1-bf1)*1e6))*2/3.14159"


;#   ifdef N15_OFFSET
;"p19=p21*atan(cnst12/abs(cnst28))*2/3.14159"
;#	else
;"p19=p21*cnst12/cnst12"
;#   endif /*N15_OFFSET*/
"p19=0.0"  ;remove flip-angle pulse 

"DELTA1=p25*2+20u"


1 ze
  d11 st0
2 3m do:f2 do:f3
3 d11
4 d11 do:f2 do:f3
   
  "p31=vd"
  "p29=d30-p31"
  "d29=d30-p31"
  
  d12 fq=cnst29:f1 fq=cnst29:f3
  d12 pl24:f1 pl0:f2 pl23:f3
  d12

;  if "d29 > 0.0"
;       {
;       2u pl24:f1
;       2u cw:f1 ph0
;       (p29 ph0):f3
;       2u do:f1
;       2u
;       }
;  else
;       {
;       8u
;       }

  if "d29 > 0.0"
       {
       4u
       (p25:sp5 ph1):f1
       4u
       2u pl24:f1 pl23:f3 fq=50000:f1 fq=50000:f3
       (center (p29 ph0):f3 (p29 ph0):f1 )
       6u fq=cnst29:f1 fq=cnst29:f3
       (p25:sp5 ph3):f1
       4u
       }
  else
       {
       DELTA1
       }
	   
  d1
  
  2u pl1:f1 pl3:f3
  d11 UNBLKGRAD
 
  (p21 ph0):f3
  4u
  p10:gp0
  4m

  (p25:sp5 ph2):f1
  4u pl1:f1
  (p1 ph0)
  2u fq=cnst30(bf ppm):f1
  (center (p11:sp1 ph11) (p11:sp2 ph12):f3)
  3u pl1:f1 pl3:f3
  (ralign (p18 ph1):f1 (p21 ph1):f3 )
  
  2u pl24:f1
  2u cw:f1 ph0
  5m 									; equilibrium delay

  if "p31 == 0.0"
       {
       if "cnst28 > 0"
	       {
	       (p19 ph5):f3
	       6u
	       (p19 ph6):f3
	       }
       else
	       {
	       (p19 ph7):f3
	       6u
	       (p19 ph8):f3
	      }
       }
  else
       {
       if "cnst28 > 0"
	      {
	       (p19 ph5):f3
	      
#   ifdef N15_OFFSET
    1u fq=cnst28:f3
#   else
    1u
#   endif /*N15_OFFSET*/
	       
	       1u pl23:f3
	       1u
	       (p31 ph0):f3
	       1u
	       1u pl3:f3
	       1u fq=cnst29:f3     
	       (p19 ph6):f3
	      }
       else
	      {
	      (p19 ph7):f3
 
#   ifdef N15_OFFSET
    1u fq=cnst28:f3
#   else
    1u
#   endif /*N15_OFFSET*/
 
	      1u pl23:f3
	      1u
	      (p31 ph0):f3
	      1u
	      1u pl3:f3
	      1u fq=cnst29:f3     
	      (p19 ph8):f3
	      }
       }
  
#   ifdef ZETA_DELAY
  (p21 ph3):f3
  d2
  (p21 ph1):f3
#   endif /*ZETA_DELAY*/ 
  
  0.1u do:f1
  2u pl1:f1
  (p18 ph3):f1
  
  2u
  p12:gp2
  d16

  (p21 ph1):f3
  2u 
  (center (p11:sp1 ph15) (p11:sp2 ph16):f3)
  3u 

  2u
  p13:gp3
  d16 fq=cnst29:f1

  3u
  (p25:sp5 ph4:r):f1
  4u pl1:f1
  (p2 ph2):f1
  4u
  (p25:sp5 ph4:r):f1
  3u
  
  p13:gp3
  d16 pl12:f2 pl16:f3
  
  2u BLKGRAD
  
  goscnp ph31 cpd2:f2 cpd3:f3

  3m do:f2 do:f3
  3m st ivd 
  lo to 2 times nbl

  3m ipp11 ipp12 ipp15 ipp16 ipp31
  lo to 3 times ns

  d11 wr #0

exit
   

ph0=0
ph1=1 
ph2=2
ph3=3

ph4=0
ph5=1
ph6=3
ph7=3
ph8=1

ph9=3 1
ph10=1 3

ph11=1 1 1 1 1 1 1 1 3 3 3 3 3 3 3 3
ph12=2 0
ph15=0 0 0 0 2 2 2 2
ph16=0 0 2 2

ph31=0 2 2 0 2 0 0 2 2 0 0 2 0 2 2 0


;***************************************

;pl0 : 120dB
;pl1 : f1 - High power pulse (default)
;pl2 : f2 - High power pulse (default)
;pl3 : f3 - High power pulse (default)
;pl12: f2 - decoupling during acqu
;pl16: f3 - decoupling during acqu
;pl23: f3 - 15N spin lockpower
;pl24: f1 - 1H spin lockpower

;sp1 : f1 - hetnuc CP [Square, ~90 Hz]
;sp2 : f3 - hetnuc CP [Square, ~90 Hz]
;sp5: f1 - 1H 90deg (Sinc.1000)

;p1 : f1 - 90 degree high power pulse
;p2 : f1 - 180 degree high power pulse
;p3 : f2 - 90 degree high power pulse
;p4 : f2 - 180 degree high power pulse
;p10: gradient - 1 ms [< 250 pts]
;p11: f1/2 - hetnuc CP [1/J(CH)]
;p12: gradient - 50 us [< 12 pts]
;p13: gradient - 500 us [< 125 pts]
;p18: f1 - atan(1Hspinlock / 1Hoffset)
;p19: f3 - atan(15Nspinlock / 15Noffset)
;p21: f3 - 90 degree high power pulse
;p22: f3 - 180 degree high power pulse
;p25: f1 - 1H 90deg (1-2 ms)

;cnst10: 15N CS diff for zeta delay [Hz]
;cnst11: 1H spinlock power [Hz]
;cnst12: 15N current spinlock power [Hz]
;cnst28: 15N offset [Hz]
;cnst29: on reson freq diff = 0 Hz
;cnst30: 1H resonance frequency [in ppm]

;d0 : incremented delay (2D)            
;d1 : recovery delay; 1-5 * T1
;d2 : zeta delay - 1/(4*cnst10)
;d11: delay for disk I/O     [30m]           
;d12: delay for power switching  [20u]
;d16: delay for gradient recovery [150u]
;p29: Tmax - Trlx
;d30: Tmax of relaxation delays
;p31: Trelax (test) - uses vdlist

;in0: 1/(2 * SW(X)) = DW(X)
;NS: 16 * n
;DS: nbl*ds >= 64
;td1: number of increments
;FnMODE: QF
;cpd2: decoupling according to sequence 
;pcpd2: f2 - 90 degree pulse for cpd
;cpd3: decoupling according to sequence 
;pcpd3: f3 - 90 degree pulse for cpd

;for z-only gradients:
;gpz0: 15%
;gpz2: 75%
;gpz3: 100%

;use gradient files:   
;gpnam0: SMSQ10.100
;gpnam2: RECT.1
;gpnam3: SMSQ10.100



                                          ;preprocessor-flags-start
;ZETA_DELAY: for suppressing signals with similar 1H but different
;			 15N resonances. zeta = 1/(4 dv), dv = cnst10 = resonance diff in hz
;			  option -DZETA_DELAY (eda:ZGOPTNS)
;N15_OFFSET: used when 15N spinlock is to be off resonance
;			  option -DC13_OFFSET (eda:ZGOPTNS)
                                          ;preprocessor-flags-end



;$Id: 15Nsel-R1rho_HC.alh, v 1.0, 19 June 2008 $

