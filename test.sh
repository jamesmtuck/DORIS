echo    "Test error analysis of fastq data..."
python3 scripts/ngs.py --range 0-26 --fastq_directory ./data/StrippedFastQ
echo    "Test one run of primer study with codeword size of 6 ..."
python3 scripts/primer_vs_density.py --codeword-size 6 --library-size=0 --samples 1000 --csv density_cw6_lib0_short_t0.csv
echo    "Evaluate prior results from primer study and produce plot ..."
python3 scripts/plot-density-tradeoff.py ./data/density_results/*
