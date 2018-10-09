#!/bin/csh

@ p0 = 115
@ p120 = $p0 + 120
@ p240 = $p0 + 240

set clists = ("1 0 0 0 0 0" "0 1 0 0 0 0" "0 0 1 0 0 0" "0 0 0 1 0 0" "0 0 0 0 1 0" "0 0 0 0 0 1")
set ZQphases = ($p0 $p120 $p240 $p0 $p240 $p120)
set DQphases = ($p0 $p240 $p120 $p0 $p120 $p240)

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1904 -dspfvs 21 -grpdly 76  \
  -xN              2048  -yN              2400  \
  -xT              1024  -yT              1200  \
  -xMODE            DQD  -yMODE        Complex  \
  -xSW        10504.202  -ySW         4226.543  \
  -xOBS         700.130  -yOBS         176.051  \
  -xCAR           0.200  -yCAR          24.380  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov


foreach i (`seq 6`)
nmrPipe -in test.fid \
| nmrPipe -fn COADD -axis Y -cList $clists[$i] -time              \
| nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 $ZQphases[$i] -p1 0.00  -verb         \
| nmrPipe -fn EXT -x1 1.5ppm -xn -1ppm -sw                                    \
   -ov -out test$i.ft1
end

addNMR -in1 test1.ft1 -in2 test2.ft1 -out test.ft1
addNMR -in1 test.ft1  -in2 test3.ft1 -out test.ft1
addNMR -in1 test.ft1  -in2 test4.ft1 -out test.ft1
addNMR -in1 test.ft1  -in2 test5.ft1 -out test.ft1
addNMR -in1 test.ft1  -in2 test6.ft1 -out test.ft1
mv test.ft1 ZQ.ft1
rm test*.ft1


foreach i (`seq 6`)
nmrPipe -in test.fid \
| nmrPipe -fn COADD -axis Y -cList $clists[$i] -time              \
| nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 $DQphases[$i] -p1 0.00  -verb         \
| nmrPipe -fn EXT -x1 1.5ppm -xn -1ppm -sw                                    \
   -ov -out test$i.ft1
end

addNMR -in1 test1.ft1 -in2 test2.ft1 -out test.ft1
addNMR -in1 test.ft1  -in2 test3.ft1 -out test.ft1
addNMR -in1 test.ft1  -in2 test4.ft1 -out test.ft1
addNMR -in1 test.ft1  -in2 test5.ft1 -out test.ft1
addNMR -in1 test.ft1  -in2 test6.ft1 -out test.ft1
mv test.ft1 DQ.ft1
rm test*.ft1


foreach i ("ZQ" "DQ")
nmrPipe -in $i.ft1 \
| nmrPipe  -fn PS -p0 0 -p1 0.00 -di -verb         \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn LP -fb \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 0.00 -p1 180.00 -di -verb         \
| nmrPipe  -fn TP                                     \
   -ov -out $i.ft2
end

