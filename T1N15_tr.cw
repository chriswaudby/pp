; 15N-T1 relaxation experiment with TROSY read-out
; for 15N, 15N13C, 2H15N and 2H15N13C labelled proteins
; written by NL 10/25/11
; see footnotes

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

;#define LABEL_CN  ; switch on for 13C labelled samples
                   ; ALSO SET 13C DECOUPLING (AT TOP, 77, AND BOTTOM OF PROGRAM)
#define TEMP_COMPENSATION


"in0=inf1*0.5"

# ifdef LABEL_CN
"d0=97u-p4*2+p7*0.66-p1*0.5"
"d25=20m-p15*0.5-p4*4-8u"
#else
"d0=100u+p7*0.66-p1*0.5"
"d25=20m-p15*0.5"
#endif /*LABEL_CN*/

"d11=30m"
"DELTA=2.65m"
"DELTA1=2.65m"
"DELTA2=2.65m-p22-p11-300u"
"DELTA3=2.65m-p23-p10-300u"
"DELTA4=260u-p24-p1*0.66"

"d27=p24+35u"

"l1=1"
"l2=1"


"cnst21=176"
"cnst22=56"
"cnst18=-800"



"spoff4=bf2*((cnst22-cnst21)/1000000)"



1       ze
        1m 
2       d11 ;do:f2                ; 13C DECOUPLING SWITCH OFF!
        1m LOCKH_OFF
        3m 
3       1m 
        1m 
4       3m 
5       2m BLKGRAD
        10u pl1:f1
        10u ;pl4:f2               ; 13C POWER
        10u pl7:f3
        (p7 ph0):f3
        10u 


;---------temperature compensation and d1 recovery delay--------- 
# ifdef TEMP_COMPENSATION

"d17=d1-p18"

        10u fq=cnst18(bf ppm):f3   
        10u pl8:f3
        (p18 ph0):f3    ; 15N pulse is applied far off-resonance      
        10u
        10u fq=0:f3
        d17
# else
        d1
# endif
        1m UNBLKGRAD
        10u pl7:f3

;------- kill steady state 15N ------------
        (p7 ph0):f3
        5u
        p20:gp6
        200u

;------- first INEPT Hz-> 2HxNz -----------
        (p1 ph0):f1
        5u
        DELTA gron0 ; soft gradient to prevent radiation damping 
        5u  groff
        (center(p1*2 ph0):f1 (p7*2 ph0):f3)
        5u
        DELTA gron0
        5u groff
 
;------- rephase  2HxNz to Nz------ --------
        (p1 ph5):f1  (p7 ph0):f3 
        5u 
        DELTA1 gron1 ; soft gradient to prevent radiation damping
        5u groff
        (center (p1*2 ph0):f1 (p7*2 ph0):f3)
        5u
        DELTA1 gron1
        5u groff
        (p7 ph6):f3 ; phase-cycle Nz, -Nz for Freeman-Hill decay
        5u
;--------------------------------------------
        (p1 ph2):f1 ; purge pulse to kill any residual HzNz
        5u
        p21:gp7  ; cleaning gradient
        100u
        100u pl0:f1
;------15N T1 relaxation period--------------

if "l2==1" goto 77  ; jump to 77 for first relaxation data point, needs to be 0 in vclist

70       d25*0.5 
#  ifdef LABEL_CN
         3u pl4:f2 
         (p4*2 ph0 3u 3u pl2 p4*2:sp4 ph0):f2
#  endif
         d25*0.5
         (p15:sp5 ph0):f1
         d25*0.5
#  ifdef LABEL_CN 
         3u pl4:f2
         (p4*2 ph0 3u 3u pl2 p4*2:sp4 ph0):f2
# endif
        d25*0.5 
        lo to 70 times c   ; delay=c*2*d25 (20ms)

;------Echo/ Anti-echo encoding for TROSY read-out------------   
 77     3u
        3u ;pl4:f2               ; 13C POWER
        3u pl1:f1 
        if "l1==1" 
        {
        (p7 ph7):f3
        10u
        p25:gp5
        200u 
        (p7*2 ph7):f3
        10u
        p25:gp5*-1
        }
        else 
        {
        (p7 ph17):f3
       	10u
        p25:gp5*-1
        200u 
        (p7*2 ph17):f3
        10u
        p25:gp5
        }
;------ t1 (15N) evolution period ------------------------------
        d0
# ifdef LABEL_CN
        (p4*2 ph0 3u 3u pl2 p4*2:sp4 ph0):f2
# endif
        d0
;------ start TROSY read-out------------------------------------
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
;goto 999 ; optimization of  water supression
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
;goto 999 ; for optimization of  water supression
        DELTA3
        p23:gp3           
        200u 
        100u pl10:f1
        (center(p10 ph10:r 5u pl1 p1*2 ph0 5u pl10 p10 ph10:r):f1 (p7*2 ph0 d27):f3)
        5u
;goto 999 ; for optimization of  water supression
        p23:gp3
        DELTA3
        DELTA4 
        (p7 ph0):f3
        5u
        p24:gp4 ; Echo/Anti-echo decoding gradient
999     5u
        5u ;pl31:f2             ; 13C DECOUPLING POWER
        20u BLKGRAMP
        go=2 ph31 ;cpds2:f2     ; 13C DECOUPLING
        1m ;do:f2               ; 13C DECOUPLING SWITCH OFF!
        1m LOCKH_OFF
        d11 wr #0 if #0 zd
        1m ivc
        1m iu2
        lo to 3 times l6
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
ph6=1 1 1 1  3 3 3 3  
ph10=2
ph11=3
ph12=0
ph13=1
ph7=1 0 3 2  
ph17=1 2 3 0 
ph31=1 2 3 0  3 0 1 2  


;-------------NOTES----------------------

;o1p = 4.7 ppm
;o2p=176 ppm (CO)
;o3p=119 ppm
 
;NS=8*n
;in0=inf/2
;SW=1/(2*in0)
;echo-antiecho in N15 (process as Complex in NmrDraw before splitting the spectra)

; loop counters
;l3: number of complex points (td1 / 2)
;l6: number of relaxation points

; 1H pulses

;p1:  90 deg hard 1H pulse @pl1
;pl1: 1H 90 deg
;pl0: 120 dB
;p10: 1000u (@ 700 MHz) 90 deg soft rectangular water flip-back pulse (pl10)
;p11: 1600u (@ 700 MHz) 90 deg Sinc1.1000 water flip-back pulse (sp11,sp12)
;p15: 1700u (@ 700 MHz) 180 deg IBurp2 pulse on 1H (sp5)
;sp5: 180 deg IBurp2 pulse on 1H (p15) 
;sp11: 90 deg Sinc1.1000 water flip-back pulse
;sp12: 90 deg Sinc1.1000 water flip-back pulse
;spnam5: IBurp2
;spnam11: Sinc1.1000
;spnam12: Sinc1.1000
;spoffs5: 2730Hz @ 700 MHz (8.6 ppm) , should be centered in amide region but not touch the water 

; 13C pulses

;p4: 13CO selective 90 deg (20.31us @ 700 MHz) @pl4
;pl2: 120 dB
;pl4: 13C 90 deg 

;sp4: 13CA selective 90 deg (p4)
;CPDPRG2: garp (aq C' decoupling)
;pcpd5: C' decoupling (140u or 280u @pl31)
;pl31: C' decoupling power

;15N pulses
;p7 : 90 deg hard 15N pulse @pl7
;p18 : maximum duration of spin-lock; temperature compensation
;pl7 :15N 90 deg
;pl8: 15N spin-lock power

; gradients
;p20: 1000u
;p21: 200u
;p22: 300u
;p23: 1000u
;p24: 60.8u Echo/Anti-echo decoding gradient
;p25: 300u Echo/Anti-echo half-encoding gradient

;for z-only gradients
;gpz0: 3%
;gpz1: 2%
;gpz2: 10%
;gpz3: 50%
;gpz4: 33%
;gpz5: -33%
;gpz6: 30%
;gpz7: -50%

;gpnam2 SINE.10
;gpnam3 SINE.50
;gpnam4 SINE.10
;gpnam5 SINE.10
;gpnam6 SINE.10
;gpnam7 SINE.10
