;19F on-resonance R1rho relaxation
;based on Overbeck (2020)
;with improved temperature compensation

/*
|
| On-resonance 19F R1rho as pseudo-3D
| with different SL lenghts read in via VPLIST
| and different SL powers read in via VALIST
|
| using hard pulses for flipdown/flipback
|
| Pseudo-3D
| Jan Overbeck
| 2020
|
*/

/*--------------------------------
; Parameters to set
; -------------------------------*/
;cnst28 : offset of SL in ppm
;p30 : maximum SL length
;p31 : heating compensation SL length
;p32 : spin lock lenght T_ex
;pl25 : spin lock power, = sp4
;VPLIST : list of spin lock lengths
;VALIST : list of spin lock powers !in dB!

#include <Avance.incl>
#include <Grad.incl>
#include <Delay.incl>

define list<pulse> plength = <$VPLIST>
define list<power> list1 = <$VALIST>

"p2=p1*2"
"d11=30m"
"l2=0"
"l3=0"
aqseq 312

1 ze
2 30m
/*--------------------------------
; calculate SL delays
; -------------------------------*/
 "p32=plength[l2]"


/* ---------------------------------
; relaxation delay (d1)
; --------------------------------*/
 d1


/* ---------------------------------
; heating compensation
; --------------------------------*/
"p31=p30-p32*pow(10,(pl30 - list1)/10)"

if "p31 > 0.0"
 {
 1u fq=cnst30(bf ppm):f1
 1u pl30:f1
 (p31 ph1):f1
 }


/* ---------------------------------
; transfer to theta and SL
; --------------------------------*/
 30m
 1u fq=cnst28(bf ppm):f1
if "p32 == 0.0"
 {
 1u pl1:f1
 p1 ph4
 }
else
 {
 1u pl1:f1
 p1 ph4
 1u list1:f1
 (p32 ph1):f1 ; <-- this is the Spin Lock
}
;-----------------------------------

/* ---------------------------------
; transfer back to z
; --------------------------------*/
 1u pl1:f1
 p1 ph5
;------------------------------------

/* ---------------------------------
; anti-ringing
; --------------------------------*/
 1u pl1:f1
 p1 ph1
 4u
 p1 ph2
 4u
 p1 ph3
;------------------------------------

; 4u BLKGRAD
 go=2 ph31
 30m mc #0 to 2
   F1QF(calclc(l2,1))
   F2QF(calclist(list1,1))
;exit
HaltAcqu, 1m
exit


ph1=0
ph2=2 0
ph3=0 0 2 2 1 1 3 3
ph4=1
ph5=3
ph31=0 2 2 0 1 3 3 1
;pl1 : f1 channel - power level for pulse (default)
;p1 : f1 channel - 90 degree high power pulse
;p2 : f1 channel - 180 degree high power pulse
;d1 : relaxation delay; 1-5 * T1
;d11: delay for disk I/O [30 msec]
;ns: 8 * n
;ds: 128

;p30: maximum SL length used for highest power (for T compensation)
;pl30: maximum SL power (for T compensation)
;cnst30: offset of SL in ppm (for T compensation) [250 ppm]




