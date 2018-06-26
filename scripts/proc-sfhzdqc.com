#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1560 -dspfvs 21 -grpdly 76  \
  -xN              3072  -yN               256  \
  -xT              1536  -yT               128  \
  -xMODE            DQD  -yMODE        Complex  \
  -xSW        12820.513  -ySW         4863.813  \
  -xOBS         800.254  -yOBS          81.098  \
  -xCAR           4.773  -yCAR         119.575  \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

nmrPipe -in test.fid \
| nmrPipe -fn COADD -axis Y -cList 1 0         \
| nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -94.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -alt                               \
| nmrPipe  -fn PS -p0 -72.00 -p1 180.00 -di -verb         \
   -ov -out test1.ft2

nmrPipe -in test.fid \
| nmrPipe -fn COADD -axis Y -cList 0 1         \
| nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -94.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -alt                               \
| nmrPipe  -fn PS -p0 18.00 -p1 180.00 -di -verb         \
   -ov -out test2.ft2

addNMR -in1 test1.ft2 -in2 test2.ft2 -add -out sum.ft2
addNMR -in1 test1.ft2 -in2 test2.ft2 -sub -out diff.ft2
