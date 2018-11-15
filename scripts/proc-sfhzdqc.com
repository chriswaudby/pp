#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1312 -dspfvs 21 -grpdly 76  \
  -xN              3072  -yN                256  \
  -xT              1536  -yT                128  \
  -xMODE            DQD  -yMODE        Complex  \
  -xSW        15243.902  -ySW         2985.075  \
  -xOBS         950.454  -yOBS          96.319  \
  -xCAR           4.771  -yCAR         116.576  \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

nmrPipe -in test.fid \
| nmrPipe  -fn MAC -all -noRd -noWr -macro shufZQ.M \
| nmrPipe  -fn SOL \
| nmrPipe -fn EM  -lb 9.0 -c 1.0                               \
#| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -zf 2                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 97.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn EXT -x1 6ppm -xn 10.5ppm -sw \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn LP -fb \
| nmrPipe  -fn EM -lb 15 -c 1.0 \
#| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -neg                               \
| nmrPipe  -fn PS -p0 -115.00 -p1 180.00 -di -verb         \
   -ov -out hzqc.ft2

nmrPipe -in test.fid \
| nmrPipe  -fn MAC -all -noRd -noWr -macro shufDQ.M \
| nmrPipe  -fn SOL \
| nmrPipe -fn EM  -lb 9.0 -c 1.0                               \
#| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -zf 2                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -83.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn EXT -x1 6ppm -xn 10.5ppm -sw \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn LP -fb \
| nmrPipe  -fn EM -lb 15 -c 1.0 \
#| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -neg                               \
| nmrPipe  -fn PS -p0 -65.00 -p1 180.00 -di -verb         \
   -ov -out hdqc.ft2

