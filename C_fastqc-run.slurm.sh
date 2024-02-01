#!/bin/bash

#   fastqc-run.slurm

#SBATCH --job-name=fastqc
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --array=1-156

echo "${SLURM_ARRAY_TASK_ID}"

###--------------- ACTUAL RUN TIME ---------------###

# Run date: 30/04/2023
# Time: 00:17:37
# Memory: 285.80 MB

###--------------- RUN AS ---------------###

# Run from terminal as
#  "sbatch fastqc-run.slurm array-fastqc" 

#-----------------ErrorCodes---------------------------#

set -o errexit
set -o pipefail
set -o nounset

#-----------------LoadModules--------------------------#

module load fastqc/0.11.9-java-11


#-----------------InputArguments-----------------------#

INPUT=(`head -n ${SLURM_ARRAY_TASK_ID} ${1} | tail -n 1`)

FILE=$(  echo ${INPUT[@]} | cut -f 2 -d " ") 

#-----------------PathVariables-----------------------#

indir=/home/mdavidson/projects/snp-chip/agrf-illumina-data
outdir=/home/mdavidson/projects/snp-chip/agrf-illumina-data/fastQC-results

#-----------------RunProgram--------------------------#

fastqc ${indir}/${FILE} -o ${outdir}


