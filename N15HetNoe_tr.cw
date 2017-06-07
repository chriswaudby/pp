; 15N HetNoe experiment with TROSY read-out
; for 15N, 2H15N, 15N13C and 2H15N13C labelled proteins
; written by NL 10/25/11
; see footnotes

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

;#define LABEL_CN  ; switch on for 13C labelled samples (NB commented out lines)


"in0=inf1*0.5"

# ifdef LABEL_CN
"d0=97u-p4*2+p7*0.66-p1*0.5"
;"d25=20m-p15*0.5-p4*4-5u"
#else
"d0=100u+p7*0.66-p1*0.5"
;"d25=20m-p15*0.5"
#endif /*LABEL_CN*/

"d11=30m"
"DELTA=2.65m-p25-200u-p10"
"DELTA1=2.65m"
"DELTA2=2.65m-p22-p11-300u"
"DELTA3=2.65m-p23-p10-310u"
"DELTA4=260u-p24-p1*0.37"

"d27=p24+35u"

"l1=1"
"l2=1"


;"cnst21=176"
;"cnst22=56"
"cnst23=8.5"

;"spoff4=bf2*((cnst22-cnst21)/1000000)"
"d8=p1"

aqseq 321

1       ze
        1m 
2       d11 ;do:f2
        1m LOCKH_OFF
        3m 
3       1m 
        1m 
4       3m 
5       2m BLKGRAD
        10u pl1:f1
        10u ;pl4:f2
        10u pl7:f3

        (p7 ph0):f3     ; 90N pulse before d1
        d1

;--------------- 1H saturation period-----------------------
;1s
;goto 10

			10u fq=cnst23(bf ppm):f1       
        
			if "l2 == 1" goto 9
8    	11m
      (p1*2 ph0):f1
      11m
			lo to 8 times l8

goto 10
 
9    	11m
      d8*2
     	11m
			lo to 9 times l8


10     10u fq=0:f1
       10u UNBLKGRAD

;------Echo/ Anti-echo encoding for TROSY read-out------------   
         if "l1==1" 
        {
        (p7 ph7):f3
        10u
        DELTA pl14:f1 
        p25:gp5
        200u
        (center(p10 ph14:r 3u pl1 p1 ph0 3u p1*2.3 ph1 3u p1 ph0 3u pl14 p10 ph15:r):f1 (p7*2 ph7):f3)
;goto 999 ;(optimisation of p10@pl14, and phase correction phcor[15])
        10u
        p25:gp5*-1
        DELTA
        }
        else 
        {
        (p7 ph17):f3
       	10u
        DELTA pl14:f1 
        p25:gp5*-1
        200u
        (center(p10 ph14:r 3u pl1 p1 ph0 3u p1*2.3 ph1 3u p1 ph0 3u pl14 p10 ph15:r):f1 (p7*2 ph17):f3)
        10u
        p25:gp5
        DELTA
        }

;------ t1 (15N) evolution period ------------------------------
        d0       
# ifdef LABEL_CN
        (p4*2 ph0 3u 3u pl2 p4*2:sp4 ph0):f2
# endif
        d0
;------ start TROSY read-out------------------------------------
        3u pl1:f1
        if "l1==1" 
        {
        (p1 ph1):f1    ; Echo 
        3u
        3u pl0:f1  
        (p11:sp11 ph11:r):f1 
        6u
        }
        else
        {
        (p1 ph3):f1   ; Anti-Echo 
        3u
        3u pl0:f1  
        (p11:sp11 ph13:r):f1 
        6u
        }
        5u pl1:f1
;goto 999 ; optimization of  water supression (p11 @ sp11)
        DELTA2
        p22:gp2
        300u
        (center (p1*2 ph0):f1 (p7*2 ph0):f3)
        7u
        p22:gp2
        DELTA2
        300u pl0:f1  
;-------------------------------------------------  
        (p11:sp12 ph12:r):f1  
        5u 
        3u pl1:f1
        if "l1==1" 
        {
        (p1 ph0):f1 (p7 ph1):f3  ; Echo
        }
        else
        {
        (p1 ph0):f1 (p7 ph3):f3     ; Anti-Echo
        }
;goto 999 ; for optimization of  water supression (p11 @ sp12)
        DELTA3
        p23:gp3           
        200u 
        100u pl10:f1
       (center(p10 ph10:r 5u pl1 p1*2 ph0 5u pl10 p10 ph10:r):f1 (p7*2 ph0 d27):f3)
        5u
        5u pl1:f1
;goto 999 ; for optimization of  water supression (p10 @ pl10)
        p23:gp3
        DELTA3
        DELTA4 
        (p7 ph0):f3
        5u
        p24:gp4 ; Echo/Anti-echo decoding gradient
999     5u
        5u ;pl31:f2
        20u BLKGRAMP
        go=2 ph31 ;cpds2:f2
        1m ;do:f2
        d11 wr #0 if #0 zd
        1m LOCKH_OFF
        2m iu2
        lo to 3 times 2
        1m iu1 
        1m igrad EA
        1m ru2
        lo to 4 times 2 
        1m id0
        1m ru1
        lo to 5 times l3
1m
1m BLKGRAD
exit    
        
ph0=0  
ph1=1  
ph2=2
ph3=3
ph5=1  
ph6=1   
ph10=2
ph11=3
ph12=0
ph13=1
ph14=2
ph15=0
ph7=1 0 3 2  
ph17=1 2 3 0 
ph31=1 2 3 0  


;-------------NOTES----------------------

;o1p = 4.7 ppm
;o2p=176 ppm (CO)
;o3p=119 ppm
 
;NS=8*n
;in0=inf/2
;SW=1/(2*in0)
;echo-antiecho in N15 (process as Complex in NmrDraw before splitting the spectra)


; 1H pulses

;p1:  90 deg hard 1H pulse @pl1
;pl1: 1H 90 deg
;pl0: 120 dB
;p10: 1000u (@ 700 MHz) 90 deg soft rectangular water flip-back pulse (pl10, pl14)
;p11: 1600u (@ 700 MHz) 90 deg Sinc1.1000 water flip-back pulse (sp11,sp12)
;sp11: 90 deg Sinc1.1000 water flip-back pulse
;sp12: 90 deg Sinc1.1000 water flip-back pulse
;spnam11: Sinc1.1000
;spnam12: Sinc1.1000

; 13C pulses

;p4: 13CO selective 180 deg (23.7*2us @ 600 MHz) @pl4
;pl4: 13C 90 deg 
;sp4: 13CA selective 180 deg (23.7*2us @ 600 MHz)  
;CPDPRG2: garp (aq C' decoupling)
;pcpd5: C' decoupling (140u or 280u @pl31)
;pl31: C' decoupling power

;15N pulses
;p7 : 90 deg hard 15N pulse @pl7
;pl7 :15N 90 deg

; gradients
;p22: 300u
;p23: 1000u
;p24: 60.8u Echo/Anti-echo decoding gradient
;p25: 300u Echo/Anti-echo half-encoding gradient

;for z-only gradients
;gpz2: 10%
;gpz3: 50%
;gpz4: 33%
;gpz5: -33%

;gpnam2 SINE.10
;gpnam3 SINE.50
;gpnam4 SINE.10
;gpnam5 SINE.10
