.PHONY: all

trial:=0
PYTHON = python3

all:
	for libsize in 0 1000000 ; do \
		for cw in 6 7 8 9 10 11 12 ; do \
			$(PYTHON) ./primer_vs_density.py --codeword-size $$cw --library-size=$$libsize --samples 10000 --samples 15000 --samples 25000 --samples 50000 --samples 100000 --samples 200000 --samples 100000 --samples 100000 --samples 100000 --samples 100000 --samples 100000 --samples 100000 --samples 200000 --samples 200000 --samples 200000 --samples 200000 --samples 200000 --csv density_cw$${cw}_lib$${libsize}_short_t${trial}.csv --verbose ; \
		done ; \
	done 

# density:
# 	for libsize in 0 1000000 ; do \
# 		for cw in 6 7 8 9 10 11 12 ; do \
# 			bsub -q "tuck serial" -J cw$${cw}lib$${libsize}_t${trial} -W 300:00 -o density_cw$${cw}lib$${libsize}_t${trial}.out ./primer_vs_density.py --codeword-size $$cw --library-size=$$libsize --samples 10000 --samples 15000 --samples 25000 --samples 50000 --samples 100000 --samples 200000 --samples 100000 --samples 100000 --samples 100000 --samples 100000 --samples 100000 --samples 100000 --samples 200000 --samples 200000 --samples 200000 --samples 200000 --samples 200000 --csv density_cw$${cw}_lib$${libsize}_short_t${trial}.csv --verbose ; \
# 		done ; \
# 	done 
