#!/bin/bash

#   fastp-array.slurm

#SBATCH --job-name=fastp
#SBATCH --mem=6G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --array=1-78

echo "${SLURM_ARRAY_TASK_ID}"

###--------------- ACTUAL RUN TIME ---------------###

# Run date: 04/05/2023
# Time: 00:00:00
# Memory: MB

###--------------- RUN AS ---------------###

# Run from terminal as
#  "sbatch fastp-array.slurm array-fastp" 

#-----------------ErrorCodes---------------------------#

set -o errexit
set -o pipefail
set -o nounset

#-----------------LoadModules--------------------------#

module load fastp/0.23.2


#-----------------InputArguments-----------------------#

INPUT=(`head -n ${SLURM_ARRAY_TASK_ID} ${1} | tail -n 1`)

R1=$(  echo ${INPUT[@]} | cut -f 2 -d " ") 
R2=$(  echo ${INPUT[@]} | cut -f 3 -d " ") 

ID=$(  echo ${INPUT[@]} | cut -f 4 -d " ") 

#-----------------PathVariables-----------------------#

indir=/home/mdavidson/projects/snp-chip/agrf-illumina-data
outdir=/home/mdavidson/projects/snp-chip/agrf-illumina-data/fastp-results

#-----------------RunProgram--------------------------#

fastp --in1 ${indir}/${R1} --in2 ${indir}/${R2} --out1 ${outdir}/${ID}_R1.trimmed.fastq.gz --out2 ${outdir}/${ID}_R2.trimmed.fastq.gz -l 50 -h ${outdir}/${ID}_wgs.html &> ${outdir}/${ID}_wgs.log


