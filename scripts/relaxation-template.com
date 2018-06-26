#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -AMX -decim 1904 -dspfvs 20 -grpdly 67.9868774414062  \
  -xN              2048  -yN                21  -zN               128  \
  -xT              1024  -yT                21  -zT                64  \
  -xMODE            DQD  -yMODE           Real  -zMODE  Echo-AntiEcho  \
  -xSW        10504.202  -ySW           21.000  -zSW         1383.892  \
  -xOBS         700.273  -yOBS           1.000  -zOBS          70.966  \
  -xCAR           4.654  -yCAR           0.000  -zCAR         115.960  \
  -xLAB              HN  -yLAB             TAU  -zLAB             15N  \
  -ndim               3  -aq2D          States                         \
| nmrPipe -fn TP  | nmrPipe -fn ZTP  | nmrPipe -fn TP -hyper \
| pipe2xyz -out ./fid/test%03d.fid -verb -ov

xyz2pipe -in fid/test%03d.fid -x -verb              \
#| nmrPipe  -fn SOL                                  \
| nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5  \
| nmrPipe  -fn ZF -auto                             \
| nmrPipe  -fn FT                                   \
| nmrPipe  -fn PS -p0 0 -p1 0.0 -di               \
#| nmrPipe  -fn EXT -xn 10ppm -x1 6ppm -sw           \
| nmrPipe  -fn TP                                   \
| nmrPipe  -fn LP -fb                               \
| nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 1  \
| nmrPipe  -fn ZF -auto                             \
| nmrPipe  -fn FT -auto                                  \
| nmrPipe  -fn PS -p0 -90.0 -p1 180.0 -di               \
| pipe2xyz -out ft/test%03d.ft2 -y

#cp ft/test001.ft2 first-plane.ft2
#sethdr first-plane.ft2 -ndim 2 -zN 1 -zT 1
#pipe2ucsf first-plane.ft2 first-plane.ucsf

