#!/bin/csh

# NB process as 2D (ymod complex)
./fid.com

# split into components
python proc-15n-ccr.py

# process each component
set expts = (naa nab nba nbb haa hab hba hbb)

foreach i (`seq 8`)

echo {$expts[$i]}.fid

nmrPipe -in {$expts[$i]}.fid \
| nmrPipe -fn SOL                                     \
#| nmrPipe -fn EM -lb 10 -c 0.5 \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 0.5    \
#| nmrPipe -fn GM -g1 8 -g2 20 -c 0.5 \
| nmrPipe  -fn ZF -zf 2                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -156.00 -p1 123.00 -di -verb       \
| nmrPipe  -fn EXT -x1 7.5ppm -xn 9ppm -sw             \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn LP -fb -ps0-0                          \
| nmrPipe -fn GM -g1 0 -g2 15 -c 0.5 \
#| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 0.5    \
| nmrPipe  -fn ZF -zf 1                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 0.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn TP                                     \
   -ov -out {$expts[$i]}.ft2

rm {$expts[$i]}.fid

end

#
# # fit intensities
# cp ../??.tab .
# autoFit.tcl -specName naa.ft2 -inTab aa.tab -outTab naa.tab
# autoFit.tcl -specName nab.ft2 -inTab ab.tab -outTab nab.tab
# autoFit.tcl -specName nba.ft2 -inTab ba.tab -outTab nba.tab
# autoFit.tcl -specName nbb.ft2 -inTab bb.tab -outTab nbb.tab
# autoFit.tcl -specName haa.ft2 -inTab aa.tab -outTab haa.tab
# autoFit.tcl -specName hab.ft2 -inTab ab.tab -outTab hab.tab
# autoFit.tcl -specName hba.ft2 -inTab ba.tab -outTab hba.tab
# autoFit.tcl -specName hbb.ft2 -inTab bb.tab -outTab hbb.tab
#
# awk '(NR>17){print $1,$18,$19}' naa.tab | sort > naa.txt
# awk '(NR>17){print $1,$18,$19}' nab.tab | sort > nab.txt
# awk '(NR>17){print $1,$18,$19}' nba.tab | sort > nba.txt
# awk '(NR>17){print $1,$18,$19}' nbb.tab | sort > nbb.txt
# awk '(NR>17){print $1,$18,$19}' haa.tab | sort > haa.txt
# awk '(NR>17){print $1,$18,$19}' hab.tab | sort > hab.txt
# awk '(NR>17){print $1,$18,$19}' hba.tab | sort > hba.txt
# awk '(NR>17){print $1,$18,$19}' hbb.tab | sort > hbb.txt
#
# join naa.txt nba.txt | sort -n | awk '{print $0, $2/$4,($2/$4)*sqrt(($3/$2)^2+($5/$4)^2)}' > N_alpha.txt
# join nab.txt nbb.txt | sort -n | awk '{print $0, $2/$4,($2/$4)*sqrt(($3/$2)^2+($5/$4)^2)}' > N_beta.txt
# join haa.txt hba.txt | sort -n | awk '{print $0, $2/$4,($2/$4)*sqrt(($3/$2)^2+($5/$4)^2)}' > H_alpha.txt
# join hab.txt hbb.txt | sort -n | awk '{print $0, $2/$4,($2/$4)*sqrt(($3/$2)^2+($5/$4)^2)}' > H_beta.txt
