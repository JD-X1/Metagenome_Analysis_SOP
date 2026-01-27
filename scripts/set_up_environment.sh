#!/bin/bash 

# Set up environment script
# this script pulls necessary Singularity images
# and downloads RQC filter data

# These are hefty downloads 120GB+ so ensure you have
# sufficient disk space and bandwidth

# Assuming you're working on


module load singularity

singularity pull docker://bryce911/bbtools:38.86
singularity pull docker://bryce911/spades:3.15.2

mkdir rqc_data; cd rqc_data; wget -O - http://portal.nersc.gov/dna/metagenome/assembly/rqcfilter/RQCFilterData.tar | tar -xf - ; cd ..
