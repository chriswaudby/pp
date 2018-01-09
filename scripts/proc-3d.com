#!/bin/csh

rm ft/*
rm *.dat
rm *.ft

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -AMX -decim 1904 -dspfvs 20 -grpdly 67.9868774414062  \
  -xN              1024  -yN                64  -zN               128  \
  -xT               512  -yT                32  -zT                64  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  -zMODE    States-TPPI  \
  -xSW        10504.202  -ySW         2270.663  -zSW        12323.829  \
  -xOBS         700.133  -yOBS          70.952  -zOBS         176.055  \
  -xCAR           4.917  -yCAR         117.214  -zCAR          41.884  \
  -xLAB              HN  -yLAB             15N  -zLAB             13C  \
  -ndim               3  -aq2D          States                         \
  -out ./ft/test%03d.fid -verb -ov


xyz2pipe -in ft/test%03d.fid -x -verb            \
| nmrPipe -fn SOL \
| nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5  \
| nmrPipe -fn ZF -auto                              \
| nmrPipe -fn FT -auto                              \
| nmrPipe -fn PS -p0 -15 -p1 0 -di -verb         \
| nmrPipe -fn EXT -x1 5.5ppm -xn 11ppm -sw              \
| pipe2xyz -out ft/test%03d.ft3 -x

# comment out SP and PS lines if doing inverse transform
xyz2pipe -in ft/test%03d.ft3 -z -verb            \
#| nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5  \
| nmrPipe -fn ZF -auto                             \
#| nmrPipe -fn PS -p0 0 -p1 0   \
| nmrPipe -fn FT -auto -di                            \
| pipe2xyz -out ft/test%03d.ft3 -z -inPlace

xyz2pipe -in ft/test%03d.ft3 -y -verb            \
| nmrPipe -fn LP -fb -ps0-0   \
| nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5   \
| nmrPipe -fn ZF -auto                              \
| nmrPipe -fn FT -auto                              \
| nmrPipe -fn PS -p0 90 -p1 0.0 -di -verb        \
| pipe2xyz -out ft/test%03d.ft3 -y -inPlace

# Hilbert inverse Fourier Transformation
xyz2pipe -in ft/test%03d.ft3 -z -verb                  \
| nmrPipe -fn HT -auto \
| nmrPipe -fn FT -inv   \
| nmrPipe -fn ZF -auto -inv \
| nmrPipe -fn LP -fb \
| nmrPipe -fn SP -off 0.5 -end 0.98 -pow 1 -c 0.5         \
| nmrPipe -fn ZF -auto                                    \
| nmrPipe -fn FT                                   \
| nmrPipe -fn PS -p0 0 -p1 0 -di -verb   \
| pipe2xyz -out ft/test%03d.ft3 -z -inPlace

#xyz2pipe -in ft/test%03d.ft3 -x -verb            \
#| nmrPipe -fn POLY -auto                            \
#| pipe2xyz -out ft/test%03d.ft3 -x -inPlace

proj3D.tcl

# merge all planes to one BIG 3D FILE!!
xyz2pipe -in ft/test%03d.ft3 -out ./final.ft

# remove the *.fid raw data
rm ft/*

