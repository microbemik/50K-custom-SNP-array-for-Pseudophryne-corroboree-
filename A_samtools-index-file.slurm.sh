#!/bin/bash

# samtools-index-file.slurm

#SBATCH --job-name=sam-index
#SBATCH --mem=50M
#SBATCH --time=00:10:00
###


###--------------- RUN AS ---------------###

# run from terminal as:
    # "sbatch samtools-index-file.slurm"

###--------------- ERROR CODES ---------------###


###--------------- LOAD MODULES ---------------#

# load sam tools 
module load gcc/10.2.0 samtools/1.12 htslib/1.12

###--------------- INPUT FILES ---------------###

###--------------- PATH VARIABLES ---------------###

###--------------- RUN PROGRAM ---------------###

# index file
samtools faidx psco-genome.fasta.gz
