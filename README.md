# DORIS

Code associated with DORIS. Reference: 

Kevin N. Lin, Kevin Volkel, James M. Tuck, Albert J. Keung. "DORIS: Dynamic DNA-based information storage." doi: https://doi.org/10.1101/836429.

[![Build Status](https://travis-ci.com/jamesmtuck/DORIS.svg?token=rCvdBqMzwWyNvxxUUbSh&branch=master)](https://travis-ci.com/jamesmtuck/DORIS)
![GitHub](https://img.shields.io/github/license/jamesmtuck/DORIS)

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Setting up the development environment](#setting-up-the-development-environment)
- [License](#license)
- [Issues](https://github.com/jamesmtuck/DORIS/issues)


## Compatibility

To use this code, most versions of python 3 should work. In particular, python 3.5 through 3.8 pass in Travis-CI. The relevant package dependences are listed in the requirements.txt file. 

## Installation 

### Use your local python environment
If you already have python 3 installed on your system, the simplest thing to do is download or checkout the code from GitHub.  Then in the DORIS directory, run the following commands:

    cd path/to/DORIS
    pip install -r requirements.txt
    pip install .

### Use Docker

Install Docker on your system. Available for most versions of Windows, Linux, and MacOS.

Download or clone this repository. Then, from a command prompt:

    cd path/to/DORIS
    docker build -t doris:1.0 .
    docker run -it doris:1.0 /bin/bash

This will bring up a virtual environment where commands can be executed. 

### Re-run our analyses

Now, you're ready to run the commands. To run our NGS analysis, use this command:

    python3 scripts/ngs.py --range 0-26 --fastq_directory ./data/StrippedFastQ

To perform a Monte Carlo simulation to design primers against a given codeword size:

    python3 scripts/primer_vs_density.py --codeword-size 6 --library-size=0 --csv density_cw6_lib0_short_t0.csv

Note, these runs can take a long time.  There is a makefile in scripts that will run all of the simulations used to produce the data in the paper.

To plot the result of our previous runs:

    python3 scripts/plot-density-tradeoff.py ./data/density_results/*
