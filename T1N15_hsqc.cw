; 15N-T1 relaxation experiment with sensitivity enhanced (Rance-Kay) HSQC read-out
; relaxation time = 40ms * vclist
; F2 = QF, F1 = Echo-AntiEcho
; for 15N, 2H15N, 15N13C and 2H15N13C labelled proteins
; written by NL 10/27/11
; see footnotes

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


;#define LABEL_CN  ; switch on for 13C labelled samples
#define TEMP_COMPENSATION

"in0=inf1*0.5"
"l3=td1/2"
"l6=td2"

"d0=10u"
"d11=30m"
"d26=p7-p1"
"DELTA=2.65m"
"DELTA1=2.65m"
"DELTA2=2.65m-200u"
"DELTA3=2.65m-200u"
"DELTA4=p24-p19*0.63"
"DELTA5=2.65m-p25-200u"
"DELTA7=2.65m"
"l1=1"
"l2=1"

#  ifdef LABEL_CN 
"DELTA6=35u+p4*4"
"d25=20m-p15*0.5-p4*4-5u"
# else
"DELTA6=30u"
"d25=20m-p15*0.5"
# endif

"cnst19=8.6"
"cnst21=176"
"cnst22=56"
"cnst18=-800"

"spoff5=(bf1*cnst19)/1000000 - o1"
"spoff4=bf2*((cnst22-cnst21)/1000000)"

"acqt0=0"


1       ze
        1m 
2       2m do:f3
#ifdef LABEL_CN
        10u do:f2
#endif
        2m LOCKH_OFF
        d11
3       3m 
4       1m 
5       1m do:f3
        1m 
#ifdef LABEL_CN
        10u do:f2
#endif
        1m BLKGRAD
        10u pl1:f1
#ifdef LABEL_CN
        10u pl4:f2
#endif


;**********  echo/anti-echo encoding *******      
        
;;---------temperature compensation and d1 recovery delay--------- 
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

;------- first INEPT Hz-> 2HxNz -------------------------
        (p1 ph0):f1
        5u
        DELTA gron0 ; soft gradient to prevent radiation damping 
        5u  groff
        (center(p1*2 ph0):f1 (p7*2 ph0):f3)
        5u
        DELTA gron0
        5u groff
 
;------- rephase  2HxNz to Nz-----------------------------
        (center(p1 ph5):f1  (p7 ph0):f3) 
        5u 
        DELTA1 gron1 ; soft gradient to prevent radiation damping
        5u groff
        (center (p1*2 ph0):f1 (p7*2 ph0):f3)
        5u
        DELTA1 gron1
        5u groff
        (p7 ph6):f3 ; phase-cycle Nz, -Nz for Freeman-Hill decay
        5u
;--------------------------------------------------------
        (p1 ph2):f1 ; purge pulse to kill any residual HzNz
;goto 999        
        5u
        p21:gp7  ; cleaning gradient
        100u
        100u pl0:f1
;------15N T1 relaxation period--------------------------

if "l2==1" goto 77

70        d25*0.5
#     ifdef LABEL_CN
         3u pl4:f2 
         (p4*2 ph0 3u 3u pl2 p4*2:sp4 ph0):f2
#     endif
         d25*0.5
         (p15:sp5 ph0):f1
         d25*0.5
#     ifdef LABEL_CN 
          3u pl4:f2
         (p4*2 ph0 3u 3u pl2 p4*2:sp4 ph0):f2
#     endif
         d25*0.5 
         lo to 70 times c   ; delay=c*2*d25 (20ms)

;---------- Echo- Antiecho encoding------------------------   

77      10u pl1:f1
#ifdef LABEL_CN
        5u pl4:f2
#else
	5u
#endif
        5u pl7:f3
	(p7 ph7):f3
;goto 999
        DELTA6;compensation for d0 15N evolution
        DELTA5 
        190u
        p25:gp5*EA
        10u             
        (center (p1*2 ph0):f1 (p7*2 ph7):f3)
        10u   
        p25:gp5*EA*-1
        190u 
        DELTA5
       				

; t1 evolution ---------------------------------------------
89      d0*0.5 gron1
	5u groff
        d0*0.5 gron1*-1
        5u groff     
#     ifdef LABEL_CN
        (center (p1*2 ph0):f1 (p4*2 ph0 3u 3u pl2 p4*2:sp4 ph0):f2)
#     else
        (p1*2 ph0):f1 
#     endif
       d0
;goto 999
;--------Rance-Kay transfer back ----------------------------
      
98     (center (p1 ph0):f1 (p7 ph8):f3)  
        DELTA2 gron2
        200u groff
       (center(p1*2 ph0):f1 (p7*2 ph0):f3)
        DELTA2 gron2
        200u groff

; ---second INEPT ----------------------------------------       
	(center (p1 ph1):f1 (p7 ph1):f3)    ;DOUBLE 90
        DELTA3 gron3
        200u groff
        (center(p1*2 ph0):f1 (p7*2 ph0):f3)
        DELTA3 gron3
        200u groff
        (d26*0.5 p19 ph0:r):f1
998     245u
        DELTA4	
	(p1*2 ph0):f1
	5u
	p24:gp4
	200u

999     5u           ; jump point for water optimisation
#ifdef LABEL_CN
        5u pl30:f2 	; 13C DECOUPLING POWER
#else
	5u
#endif
        10u pl31:f3
        20u BLKGRAMP

	; acquisition
#ifdef LABEL_CN
        go=2 ph31 cpds3:f3 cpd2:f2	
        500u do:f3
        500u do:f2
#else
        go=2 ph31 cpds3:f3 
        1m do:f3
#endif
        1m LOCKH_OFF
        d11 wr #0 if #0 zd

	; loop over relaxation times
        1m ivc
        1m iu2
        lo to 3 times l6

	; loop over 15N evolution
        1m ip8*2
        1m igrad EA
	1m ru2
	lo to 4 times 2 
        1m id0
        lo to 5 times l3


1m do:f3
1m BLKGRAD
exit    
        
ph0=0
ph1=1          ;check right phase for Boltzmann !!!!!
ph2=2
ph5=1 1 1 1 
ph6=1 1 3 3
ph7= 1 3 
ph8= 0  
ph31=1 3 3 1 


;-------------NOTES----------------------

;o1p = 4.7 ppm
;o2p=176 ppm (CO)
;o3p=119 ppm
 
;NS=4*n
;in0=inf/2
;SW=1/(2*in0)
;echo-antiecho in N15 (process as Complex in NmrDraw before splitting the spectra)


; 1H pulses

;p1:  90 deg hard 1H pulse @pl1
;p19 last 90 deg hard 1H pulse @pl1, can be adjusted for improving water-supression
;pl1: 1H 90 deg 
;pl0: 120 dB
;p15: 1618us (@ 700 MHz) 180 deg IBurp2 pulse on 1H (sp5)
;sp5: 180 deg IBurp2 pulse on 1H (p15) 
;spnam5: IBurp2
;spoffs5: 8.6 ppm, should be centered in amide region but not touch the water 

; 13C pulses

;p4: 13CO selective 180 deg (23.7*2us @ 600 MHz) @pl4
;pl4: 13C 90 deg 
;sp4: 13CA selective 180 deg (23.7*2us @ 600 MHz)  


;15N pulses
;p7 : 90 deg hard 15N pulse @pl7
;d17: calculated spin-lock delay
;p18: maximum duration of spin-lock; temperature compensation
;pl7 :15N 90 deg
;pl8: 15N spin-lock power
;CPDPRG3: garp (aq 15N decoupling)
;pcpd3: 15N decoupling (200u @pl31)
;pl31: 15N decoupling power



; gradients
;p20: 1000u
;p21: 200u
;p24: 201u Echo/Anti-echo decoding gradient
;p25: 1000u Echo/Anti-echo half-encoding gradient

;for z-only gradients
;gpz10: -1.5%
;gpz2: -5%
;gpz3: -0.5%
;gpz4: 40%
;gpz5: 40%
;gpz0: 3%
;gpz1: 2%
;gpz6: 30%
;gpz7: -50%

;gpnam4 SINE.10
;gpnam5 SINE.10
;gpnam6 SINE.50
;gpnam7 SINE.10


