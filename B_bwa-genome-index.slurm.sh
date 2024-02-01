#!/bin/bash

#  bwa-genome-index.slurm

#SBATCH --job-name=bwa
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00        


###--------------- RUN AS ---------------###

# Run from terminal as
    # "sbatch bwa-genome-index.slurm" 

###--------------- ERROR CODES ---------------###

set -o errexit
set -o pipefail
set -o nounset

###--------------- LOAD MODULES ---------------###

module load gcc/10.2.0 samtools/1.12 htslib/1.12
module load bwa/0.7.17

###--------------- INPUT FILES ---------------###



###--------------- PATH VARIABLES ---------------###


###--------------- RUN PROGRAM ---------------###

bwa index psco-genome.fasta.gz




