#include <Avance.incl>
#include <Grad.incl>

define delay TIME_T2
"TIME_T2 = 40m"

define delay DEL
define delay DEL_1
define delay DEL_2

define delay TAUB
"TAUB = 960u - p24 - 2u"
define delay COMPH
"COMPH = p1"

define loopcounter COUNTER

;define list<loopcounter> ncyc = {0}
;define list<loopcounter> ncyc = {0 1 20}
;define list<loopcounter> ncyc = {0 1 40}
;define list<loopcounter> ncyc = {0 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 24 28 32 36 40}
;define list<loopcounter> ncyc = {0 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20}
;define list<loopcounter> ncyc = {0 20 1 18 2 16 3 14 4 12 5 10 6 9 7 8}
define list<loopcounter> ncyc = {0 40 1 36 2 32 3 28 4 24 5 20 6 18 7 16 8 14 9 12 10}
;define loopcounter ncyc
;"ncyc = l0" ;"ncyc=0"

;"d3=3.75m - p23 - 2u"
;"d4=3.75m - p23 - 2u - p8"
"d3=4m - p23 - 2u + 0.6366*p1"
"d4=4m - p23 - 2u - p8"
"acqt0=0"

"d11=1m"
"d12=1m"

;carbon evolution time
"in0=inf2/2"
"d0=in0*0.5 - p8*0.5*0.63662 - p1*2 - 2u - p1"

aqseq 312

1 ze 
  10m st0
2 2u do:f2
3 2u do:f2
  1m BLKGRAD
4 d1 do:f2

if "ncyc == 0"
  {
   "DEL = TIME_T2/4 - p8"
   "DEL_1 = DEL - p8*0.63662*0.5"
   "DEL_2 = DEL - p1"
   "COUNTER = 0"
  }
else
  {
   "DEL = TIME_T2/(4*ncyc) - p8"
   "DEL_1 = TIME_T2/(4*ncyc) - p8 - p8*0.63662*0.5"
   "DEL_2 = DEL - p1"
   "COUNTER = ncyc - 2"
  }

  1m UNBLKGRAD
  10u pl8:f2
  10u pl0:f1
  (p2:sp2 ph0:r):f1    ; SOLVENT SUPPRESSION FLIP-BACKS
  2u pl0:f1
  (p8 ph0):f2
  2u
  p20:gp0
  2.0m pl1:f1
  

   (p1 ph0):f1
   2u
   p23:gp3
   d3

   (p8 ph1):f2
   2u
   p24:gp4
   TAUB
   (center (p8*2 ph0):f2 (p1*2 ph0):f1)
   2u
   p24:gp4
   TAUB
   (p8 ph2):f2
;******	CPMG STARTS HERE **************
if "ncyc == 1"
{  
   DEL
   (p8*2 ph5):f2
   DEL_2
}

if "ncyc == 2"
{  
   DEL
   (p8*2 ph5):f2
   DEL
   DEL
   (p8*2 ph5):f2
   DEL_2 
}

if "ncyc > 2"
{  
   DEL
   (p8*2 ph5):f2
   DEL
10 DEL
   (p8*2 ph5):f2
   DEL 
lo to 10 times COUNTER
   DEL
   (p8*2 ph5):f2
   DEL_2 
}

  (p1*2 ph3):f1

if "ncyc == 1"
{  
   DEL_2
   (p8*2 ph5):f2
   DEL
}

if "ncyc == 2"
{  
   DEL_2
   (p8*2 ph5):f2
   DEL
   DEL
   (p8*2 ph5):f2
   DEL
}

if "ncyc > 2"
{  
   DEL_2
   (p8*2 ph5):f2
   DEL
20 DEL
   (p8*2 ph5):f2
   DEL
lo to 20 times COUNTER
   DEL
   (p8*2 ph5):f2
   DEL 
}  
;****** CPMG ENDS HERE ****************
;******	carbon evolution **************
if "ncyc == 0"
{  
   d0
   (p1 ph0 2u p1*2 ph7 2u p1 ph0):f1
   d0
}
else
{
   COMPH
   d0
   (p1 ph0 2u p1*2 ph7 2u p1 ph0):f1
   d0
   COMPH
}
;****** back to protons ***************
   (p8 ph4):f2
   2u
   p23:gp3
   d4 pl31:f2

;start looping experiment
   goscnp ph31 cpd2:f2

   3m do:f2
   3m st ncyc.inc   ; nbl = number of CPMG points
   lo to 2 times nbl

   3m ipp1 ipp2 ipp3 ipp31
   lo to 3 times ns  ; CPMG innermost loop, then ns

   1m BLKGRAD
   d1 mc #0 to 4
      F1QF()
      F2PH(ip4, id0)

;   d11 wr #0 if #0 zd
;
;   d12 ncyc.inc
;   lo to 3 times l2    ; l2 = number of CPMG points
;
;   d12 ip4
;   lo to 4 times 2     ; states-TPPI
;
;   d12 ip31
;   d12 ip31    
;   d12 id0
;   lo to 5 times l3    ; l3 = number of complex points

1m BLKGRAD
exit 

ph0=0 
ph1=0 2
ph2=1 1 3 3
ph3=0 0 0 0 1 1 1 1
ph4=0
ph5=1
ph7=1
ph31=0 2 0 2 2 0 2 0


;pl1 : 1H hard 90
;pl8 : 13C hard 90 and CPMG (19 kHz)
;pl31 : 13C decoupling (2 kHz)
;p1 : 1H hard 90
;p2 : 1000 us, water flip-down
;p8 : 13C hard 90 (13.16 us for 19 kHz)
;p20 : purge gradient [1 ms]
;p23 : coherence selection [600 us]
;p24 : coherence selection [200 us]
;sp2 : 1H water flip-down (sinc1.1000)
;spoff2 : place on residual HDO resonance
;cpd2 : 13C decoupling, WALTZ-16
;pcpd2 : 13C decoupling (2 kHz)
;gpz0 : 10%
;gpz3 : 40%
;gpz4 : -40%
;o1p : 1H carrier on methyl resonances

