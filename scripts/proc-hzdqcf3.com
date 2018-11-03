#!/bin/csh

# template for processing hzdqcf3.cw experiment

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1904 -dspfvs 21 -grpdly 76  \
  -xN              4096  -yN               1400  \
  -xT              2048  -yT               700  \
  -xMODE            DQD  -yMODE        Complex  \
  -xSW        10504.202  -ySW         2483.855  \
  -xOBS         700.133  -yOBS          70.952  \
  -xCAR           4.916  -yCAR         117.225  \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

nmrPipe -in test.fid \
| nmrPipe -fn COADD -axis Y -cList 1.0 1.0 -1.0 -1.0 -time                   \
| nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -50.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn EXT -x1 7ppm -xn 9ppm -sw \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -90.00 -p1 180.00 -di -verb         \
   -ov -out hdqc.ft2


nmrPipe -in test.fid \
| nmrPipe -fn COADD -axis Y -cList 1.0 1.0 1.0 1.0 -time                   \
| nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -50.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn EXT -x1 7ppm -xn 9ppm -sw \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -auto                               \
| nmrPipe  -fn FT -neg                                \
| nmrPipe  -fn PS -p0 90.00 -p1 180.00 -di -verb         \
   -ov -out hzqc.ft2
