;ipap-ct-hsqc.gk
;Georg Kontaxis and Ad Bax
;Submitted J. Biomol. NMR, Feb. 2001.
;adapted CW Jan 2020

prosol relations=<triple>

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>


;bodenhausen-ruben  constant time in t1
;constant time tuned to 1/Jcc
;put 13C carrier at 20 ppm, use high power pl2/pl4,
;and use hsec.3 at 43 ppm 

;pulses
;p1 = 90 high power 1H @pl1
;p10 = 90 selective 1H for WATERGATE using sinc1.0 (1ms at 800 MHz)
;p3 = 90 high power 13C @pl2
;p4*2 = 180 low power 13C null at CO (p4=17.8u 800 MHz) @pl4
;p14 = 13C inversion pulse 
;p31 = low power dec on 15N for WALTZ16 dec (p90=200u)
;p29 = low power dec on 13C for GARP dec during acquisition (p90=90u)

;power levels
;pl1 1H high power
;pl0 120dB to preset shaped pulses on 1H channel
;pl2 13C high power
;pl3 15N high power
;pl4 13C power level for semi selective 180
;pl6 120dB to preset shaped pulses on 13C channel
;pl16 15N power for cpd
;pl12 13C power for cpd

;timing compensation delays
"d26=p3-p1"
"d25=0.5*d26"
"d27=0.5*p14-p1"
"d28=p1*2.15+3u-p4"

;disk I/O
"d11=50m"
"d12=0.5m"

;INEPT delays
"d2=1.8m-p21-100u"       ;forward INEPT     
"d4=1.8m-p10-p22-100u"    ;reverse INEPT

;incremented delays for 13C t1 evolution
;set for 0 linear phase correction in F1
"d14=13.3m"
"d20=10u"
"d21=10u+10u+d14+p3*1.27"

;l3 number of t1 increments
;increments for evolution delays
;in21=d21/l3, in20=in21
;13C sw = 1/(in20+in21)
"in20=inf1/2"
"in21=in20"

;JCH evolution 'IPAP' delay
"d22=2u"

;increment for 'IPAP' delay
;in22=666u (0.5/(6J)) assuming JCH=125 Hz for methyl groups
"in22=666u"
;2*d22 will be set to 0,1/6J,1/3J,1/2J in turn for each t1 incremnt

;gradients
;gradient Pulses cannot be hardwired, set them by hand
;all gradients are applied with a sinc shape
;"p18=4.0m"  ; 30% xy
;"p19=3.5m"  ; 50% x
;"p20=3.0m"  ; 50% y
;"p21=1.0m"  ;-20% z
;"p22=0.5m"  ; 30% z
;"p28=3.0m"  ; 40% z

;compiler flags
;CARBON: if undefined, only the first t1 increment is acquired: for calibration
;IPAP:   if undefined, only regular 1H coupled 13C CT-HSQC in acquired

#define CARBON          
#define IPAP

;------------start actual sequence code ---------------

1    ze 
     1m
2    d11 do:f2 do:f3
3    d12*3
4    d12*3

# ifdef OFFRES_PRESAT
  30u fq=cnst19(bf hz):f1
  4u pl9:f1
  d1 cw:f1 ph1
  4u do:f1
  30u fq=0:f1
# else
  d1
# endif /*OFFRES_PRESAT*/

     10u pl1:f1
     10u pl2:f2
     10u pl16:f3
     10u UNBLKGRAD

;---------------- start sequence ------------------------
     (p3 ph5):f2
     p18:gp18           ;saturate 13C magnetization
     4m
;---------------- INEPT H(aliph) --> C(aliph)------------
     (p1 ph0):f1
     100u
     p21:gp21
     d2
     (d26 p1*2 ph0):f1 (p3*2 ph1):f2
     100u
     p21:gp21
     d2  
     (p1 ph11):f1
     2u
     p20:gp20           ; strong gradient, purge 2HzCz state
     2m
;----------------- end H to C INEPT ---------------------
;----------------- now transverse on C ------------------
     (p3 ph12):f2
;----------------- IPAP selection evolve JCH during d22*2
     d22 pl4:f2   ;         ;13C 180/1H 180(composite)
     (p1 ph0 3u p1*2.3 ph1 3u p1 ph0):f1 (d28 p4*2 ph0):f2
     d22 pl2:f2
;---------------- z-filter -----------------------------
     (p3 ph13):f2
     100u
     p28:gp28              ; strong gradient for z-filter
     300u
     (p3 ph2):f2
;----------------- start CT evolution on carbon (t1) -----
     5u ;cpd3:f3             ; cpd 15N during t1
     d20 ;pl6:C2
     ;(p14:sp5 ph0):f2      ; C' refocussing
     d14
     5u ;pl6:f2
     ;(p3*2 ph0):f2  ; hard 13C 180
     (p14:sp3 ph0):f2      ; C(aliph) refocussing hsec.3
     5u                    ; center t1 evolution
     d21
     5u ;pl6:f2
     (p14:sp3 ph0):f2      ; compensation for hsec.3 pulse
     5u pl2:f2
     5u ;do:f3
;----------------- end CT evolution on carbon -----------
;----------------- INEPT C --> H ------------------------
     (p3 ph3):f2
     p19:gp19              ; strong gradient, purge 2HzCz state
     (d25 p1 ph0):f1 
     100u
     p22:gp22              ; WG gradient
     d4 pl0:f1
     (p10:sp0 ph14):f1       ; WG soft 1H 90
     2u pl1:f1
     (d26 p1*2 ph4):f1 (p3*2 ph5):f2    ;water gate for H2O suppression
     2u pl0:f1
     (p10:sp0 ph14):f1       ; WG soft 1H 90
     100u 
     p22:gp22              ; WG gradient
     d4 pl12:f2 BLKGRAD
     go=2 ph31 cpd2:f2
     5u do:f2
     d11 wr #0 if #0 zd
#ifdef CARBON
     1u ip2
     lo to 3 times 2       ;STATES-TPPI
#endif
#ifdef IPAP
0.1m id22                  ;increment IPAP delay by 0.5/6J
0.1m ip13                  ;increment phase of IPAP pulse
    lo to 3 times 4        ;need 4 IPAP experiments for CH3 sepapration
#endif
#ifdef CARBON
     1u ip31
     1u ip31
d12 rd22                   ;reset IPAP delay
d12 id20                   ;CT 13C evolution
d12 dd21                   ;CT 13C evolution
lo to 4 times l3           ;l3 .. number of t1 increments
#endif
d12 do:f3
d12 do:f2
exit 
  
ph0=0 
ph1=1
ph11=1 1 1 1 3 3 3 3
ph12=0
ph13=0
ph2=0 2 
ph3=0 0 2 2
ph4=0  
ph14=2
ph5=0 
ph6=0
ph31=0 2 2 0 2 0 0 2

