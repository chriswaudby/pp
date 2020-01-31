#!/bin/csh

rm *.fid *.ft2

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1312 -dspfvs 21 -grpdly 76  \
  -xN              4096  -yN               1024  \
  -xT              2048  -yT                512  \
  -xMODE            DQD  -yMODE        Complex  \
  -xSW        15243.902  -ySW         7407.407  \
  -xOBS         950.454  -yOBS         238.996  \
  -xCAR           4.771  -yCAR          22.740  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

# splits the original ser file (test.fid) of a ch13 ipap-ct-hsqc
# with four interleaved experimients
# into four separate data sets 
# A(0), B(1/6J), C(1/3J), D(1/2J)
# Georg Kontaxis and Ad Bax, 22.Feb.2001

set matA = ( 1 0 \
             0 1 \
             0 0 \
             0 0 \
             0 0 \
             0 0 \
             0 0 \
             0 0 )

set matB = ( 0 0 \
             0 0 \
             1 0 \
             0 1 \
             0 0 \
             0 0 \
             0 0 \
             0 0 )

set matC = ( 0 0 \
             0 0 \
             0 0 \
             0 0 \
             1 0 \
             0 1 \
             0 0 \
             0 0 )

set matD = ( 0 0 \
             0 0 \
             0 0 \
             0 0 \
             0 0 \
             0 0 \
             1 0 \
             0 1 )

nmrPipe -in test.fid \
| nmrPipe -fn QMIX -ic 8 -oc 2 -cList $matA -time \
  -out ./A.fid -verb -ov

nmrPipe -in test.fid \
| nmrPipe -fn QMIX -ic 8 -oc 2 -cList $matB -time \
  -out ./B.fid -verb -ov

nmrPipe -in test.fid \
| nmrPipe -fn QMIX -ic 8 -oc 2 -cList $matC -time \
  -out ./C.fid -verb -ov

nmrPipe -in test.fid \
| nmrPipe -fn QMIX -ic 8 -oc 2 -cList $matD -time \
  -out ./D.fid -verb -ov

# lc.com
# Georg Kontaxis and Ad Bax, 26.Feb.2001
# script to generate linear combinations of ipap-ct-hsqc datasets
# for separation of ch3 multiplet into four subspectra
# note that the coefficients may need further adjustment by exp(2D/T2) 
# to account for 13C T2 relaxation



addNMR -in1 A.fid -in2 D.fid -out A+D.fid -c1 1.0 -c2 1.14 -verb
addNMR -in1 A.fid -in2 D.fid -out A-D.fid -c1 1.0 -c2 1.14 -sub -verb

addNMR -in1 B.fid -in2 C.fid -out B+C.fid -c1 1.02 -c2 1.07 -verb
addNMR -in1 B.fid -in2 C.fid -out B-C.fid -c1 1.02 -c2 1.07 -sub -verb

# most downfield line
addNMR -in1 A+D.fid -in2 B+C.fid -out A+2B+2C+D.fid -c1 1.0 -c2 2.0 -verb

# center upfield line
addNMR -in1 A+D.fid -in2 B+C.fid -out A-B-C+D.fid -sub -c1 1.0 -c2 1.0 -verb

# center downfield line
addNMR -in1 A-D.fid -in2 B-C.fid -out A+B-C-D.fid -c1 1.0 -c2 1.0 -verb

# most upfield line
addNMR -in1 A-D.fid -in2 B-C.fid -out A-2B+2C-D.fid -sub -c1 1.0 -c2 2.0 -verb

#!/bin/csh

foreach spec (A-2B+2C-D  A+2B+2C+D  A-B-C+D  A+B-C-D)
   echo Processing Spectrum $spec

   nmrPipe -in $spec.fid \
   #| nmrPipe  -fn POLY -time                           \
   | nmrPipe  -fn SOL \
   | nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 \
   | nmrPipe  -fn ZF -zf 2                             \
   | nmrPipe  -fn FT -verb                             \
   | nmrPipe  -fn PS -p0 91.0 -p1 0.0 -di             \
   | nmrPipe  -fn BASE -nw 100 -nl 8ppm -2ppm \
   | nmrPipe  -fn EXT -x1 2.2ppm -xn -0.5ppm -sw                       \
   | nmrPipe  -fn TP                                   \
   | nmrPipe  -fn LP -ps0-0 -fb \
   | nmrPipe  -fn SP -off 0.5 -end 0.98 -pow 2 -c 0.5 \
   | nmrPipe  -fn ZF -zf 2                         \
   | nmrPipe  -fn FT -alt -verb                        \
   | nmrPipe  -fn PS -p0 0.0 -p1 0.0 -di               \
#   | nmrPipe  -fn TP                                   \
#   | nmrPipe  -fn POLY -auto                           \
      -ov -out $spec.ft2
end

rm *.fid
