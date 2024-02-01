#!/bin/bash

#  <multiqc.slurm>

#SBATCH --job-name=multiqc
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=01:30:00


###--------------- RUN AS ---------------###

# Run from terminal as
    # "sbatch multiqc.slurm"

###--------------- ERROR CODES ---------------###

set -o errexit
set -o pipefail
set -o nounset

###--------------- LOAD MODULES ---------------###

module load multiqc/1.9-python-3.7.4

###--------------- INPUT FILES ---------------###



###--------------- PATH VARIABLES ---------------###

indir=/home/mdavidson/projects/snp-chip/agrf-illumina-data/fastQC-results/trimmed-files
outdir=/home/mdavidson/projects/snp-chip/agrf-illumina-data/fastQC-results/trimmed-files

###--------------- RUN PROGRAM ---------------###

multiqc /home/mdavidson/projects/snp-chip/agrf-illumina-data/fastQC-results/#!/bin/bash

#  <multiqc.slurm>

#SBATCH --job-name=multiqc
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=01:30:00


###--------------- RUN AS ---------------###

# Run from terminal as
    # "sbatch multiqc.slurm"

###--------------- ERROR CODES ---------------###

set -o errexit
set -o pipefail
set -o nounset

###--------------- LOAD MODULES ---------------###

module load multiqc/1.9-python-3.7.4

###--------------- INPUT FILES ---------------###



###--------------- PATH VARIABLES ---------------###

indir=/home/mdavidson/projects/snp-chip/agrf-illumina-data/fastQC-results/trimmed-files
outdir=/home/mdavidson/projects/snp-chip/agrf-illumina-data/fastQC-results/trimmed-files

###--------------- RUN PROGRAM ---------------###

multiqc /home/mdavidson/projects/snp-chip/agrf-illumina-data/fastQC-results/trimmed-files



