# DORIS

[![Build Status](https://travis-ci.com/jamesmtuck/DORIS.svg?token=rCvdBqMzwWyNvxxUUbSh&branch=master)](https://travis-ci.com/jamesmtuck/DORIS)
![GitHub](https://img.shields.io/github/license/jamesmtuck/DORIS)

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [License](#license)
- [Issues](https://github.com/jamesmtuck/DORIS/issues)

# Overview

Code associated with DORIS. Reference: 

Kevin N. Lin, Kevin Volkel, James M. Tuck, Albert J. Keung. *DORIS: Dynamic DNA-based information storage.* doi: https://doi.org/10.1101/836429.

# Documentation

As documentation for the softwarwe becomes available, it will be placed under the docs folder.

# System Requirements

## Hardware Requirements
DORIS requires only a standard computer with enough RAM and compute power to support the needed operations.

## Software Requirements
### OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

+ macOS: Catalina 10.15.3
+ Linux: Ubuntu 18.04.3

Note that most OSes will support our software by using Docker.

### Python Dependences

Our code has been tested on python versions 3.5 to 3.8. It has the following dependences:

```
nose
sphinx
editdistance
statistics
biopython
matplotlib
pandas
numpy
scipy
primer3-py
python-Levenshtein
```

# Installation Guide

## Use your local python environment
If you already have python 3 installed on your system, the simplest thing to do is download or checkout the code from GitHub.  Then in the DORIS directory, run the following commands:

    git clone https://github.com/jamesmtuck/DORIS
    cd DORIS
    pip install -r requirements.txt
    pip install .

## Use Docker

If you do not already have Docker, you will need to install Docker on your system. It is available for free for most versions of Windows, Linux, and MacOS. You may need to be the owner or administrator of the system to install Docker.

Instructions for setting up Docker.  From a command prompt, run these commands:

    git clone https://github.com/jamesmtuck/DORIS
    cd DORIS
    docker build -t doris:1.0 .
    docker run -it -v `pwd`:/DORIS doris:1.0 /bin/bash

This will bring up a command prompt in a Linux container where commands can be executed. 

## Run Our Analyses

Now, you're ready to run the commands. To run our NGS analysis, use this command:

    python3 scripts/ngs.py --range 0-26 --fastq_directory ./data/StrippedFastQ
    
This command should run to completion within a few minutes. The result of this script is a directory named DORIS_DATA that contains files with analysis. --range specifies the range of fastq files from the directory specified by --fastq_directory to run the analysis that the analysis will run on in one process. The order of fastq files is based on the pyton sort() order of the files in the fastq directory. After running analysis over a set of fastq files, two output files will be created for each group of files in the DORIS_DATA directory. One being *_count.csv, a csv file that contains the number of occurrences for each barcode for each analyzed fastq. A total number of reads is also provided for each fastq file. The second file is *_error_rate.csv that shows the base-position error rate for the payload region for each of the 'NNN' and 'NNNNN' barcodes for each fastq file. Position '0' of the error rates corresponds to the first base of the payload region. The prefix of each of these file types will follow that of the name of the first fastq file indicated by the --range argument.

To perform a Monte Carlo simulation to design primers against a given codeword size:

    python3 scripts/primer_vs_density.py --codeword-size 6 --library-size=0 --csv density_cw6_lib0_short_t0.csv

This command produces a CSV file named density_cw6_lib0_short_t0.csv, as specified in the command line. Note, these runs can take a long time (many hours, depending on the speed of the CPU).  

There is a makefile in scripts that will run all of the simulations used to produce the data in the paper.  To collect all data from scratch, this will take multiple days of execution time on a single CPU. We recommend running all of the jobs in parallel on an HPC cluster. 

For convenience, we are provided result from our own simulations in data/density_results. To plot the result of our previous runs:

    python3 scripts/plot-density-tradeoff.py ./data/density_results/*
    
This command runs in a few seconds and produces a file named CodewordVsDensity.png.    
   
# License

This software is released under the MIT License.


