#!/usr/bin/env julia

# processing for HQQC experiment:
#
# bruk2pipe -verb -in ./ser \
#   -bad 0.0 -ext -aswap -AMX -decim 1056 -dspfvs 21 -grpdly 76  \
#   -xN              4096  -yN               2688  \
#   -xT              2048  -yT               1344  \
#   -xMODE            DQD  -yMODE        Complex  \
#   -xSW        18939.394  -ySW         1075.500  \
#   -xOBS         950.450  -yOBS         238.995  \
#   -xCAR           0.500  -yCAR          16.950  \
#   -xLAB              1H  -yLAB             13C  \
#   -ndim               2  -aq2D         Complex  \
#   -out ./test.fid -ov
#
# ./proc.jl
#
# sethdr DQQ.fid -yN 256 -yT 128
# sethdr QQ.fid -yN 256 -yT 128
#
# foreach i ("QQ" "DQQ")
# nmrPipe -in $i.fid \
# | nmrPipe -fn MAC -macro $NMRTXT/bruk_ranceY.M -noRd -noWr     \
# | nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
# | nmrPipe  -fn ZF -auto                               \
# | nmrPipe  -fn FT -auto                               \
# | nmrPipe  -fn PS -p0 74.00 -p1 0.00 -di -verb         \
# | nmrPipe -fn EXT -x1 1.5ppm -xn -1ppm -sw                                    \
# | nmrPipe  -fn TP                                     \
# | nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 1.0    \
# | nmrPipe  -fn ZF -auto                               \
# | nmrPipe  -fn FT -auto                                \
# | nmrPipe  -fn PS -p0 -101.00 -p1 194.00 -di -verb         \
#    -ov -out $i.ft2
# end


function proc(outputname, deltap1, deltap2)
	filename = "test.fid";
	td = 2048

	# input phase cycle
	#phi1 = [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6] * 2pi / 7
	#phi2 = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2] * 2pi / 3
	phi1 = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6] * 2pi / 7
	phi2 = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2] * 2pi / 3

	nphase = 21

	# sum up echo and anti-echo pathways
	phirx1 = -deltap1*phi1 + deltap2*phi2
	phirx2 =  deltap1*phi1 - deltap2*phi2

	phirx1 = exp.(1im * phirx1)
	phirx2 = exp.(1im * phirx2)


	@show npoints = Int(filesize(filename)/4 - 512)
	@show ncomplex = Int(npoints / (2td*nphase))
	# preallocate data (and dummy header)
	header = zeros(Float32, 512)
	y = zeros(Float32, npoints)

	# read the file
	open(filename) do f
	    read!(f, header)
	    read!(f, y)
	end

	y = reshape(y, td, 2, :)
	yc = y[:,1,:] + 1im * y[:,2,:]

	y = reshape(yc, td, nphase, :)
	phirx1 = reshape(phirx1, 1, nphase, 1)
	phirx2 = reshape(phirx2, 1, nphase, 1)
	y1 = sum(y .* phirx1, dims=2)
	y2 = sum(y .* phirx2, dims=2)
	
	y1 = reshape(y1, td, ncomplex)
	y2 = reshape(y2, td, ncomplex)

	# read the file
	open(outputname, "w") do f
	    write(f, header)
	    for i=1:ncomplex
			write(f, Float32.(real.(y1[:,i])))
			write(f, Float32.(imag.(y1[:,i])))
			write(f, Float32.(real.(y2[:,i])))
			write(f, Float32.(imag.(y2[:,i])))
		end
	end

	#run(`sethdr $outputname -yN $nrelax -yT $nrelax`)
end

proc("QQ.fid", 3, 1)
proc("DQQ.fid", 3, -1)
