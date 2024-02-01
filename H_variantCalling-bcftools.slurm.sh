#!/bin/bash
#
#                             Bcftools mpileup
#
#SBATCH --job-name=mpileup
#SBATCH --mem=30G        #smaller test file used 1.23GB
#SBATCH --time=1-00:00:00  #smaller test file took 12:49:13
#SBATCH --array=2-13,15-20


#Load in terminal as "sbatch variantCalling-bcftools.slurm bcftools-array"

#----------------------Error_codes------------------#

set -o errexit
set -o pipefail
set -o nounset

#----------------------Load the module/software Environment------------------#

module load gcc/10.3.0 bcftools/1.15

#----------------------input_arguments------------------#

echo "${SLURM_ARRAY_TASK_ID}"

INPUT=(`head -n ${SLURM_ARRAY_TASK_ID} ${1} | tail -n 1`) 

ID=$(  echo ${INPUT[@]} | cut -f 2 -d " ")  #e.g. 08UCPB10

REF=/data/projects/punim1823/genome/psco-genome.fasta.gz  #reference genome with path


#----------------------Path_variables------------------#

indir=/home/mdavidson/projects/snp-chip/agrf-illumina-data/bams
outdir=/home/mdavidson/projects/snp-chip/agrf-illumina-data/vcfs

#----------------------Run_program------------------#

bcftools mpileup -f ${REF} ${indir}/${ID}.merged.nd.bam | bcftools call -mv -Oz -o ${outdir}/${ID}.raw.bcftools.vcf.gz