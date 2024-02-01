#!/bin/bash

#      1-runBWA.slurm

#SBATCH --job-name=BWA_alignment
#SBATCH --mem=30G           #PSCO tests used 24 GB
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16 # BWA appears to thread well so feel free to go up to 16/32 here.
#SBATCH --time=1-00:00:00   #PSCO tests took 16h
#SBATCH --array=1-78

# Run from terminal as "sbatch runBWA.slurm bwa-array"  


#--------------------ERROR-CODES----------------------------#

set -o errexit
set -o pipefail
set -o nounset

#--------------------ENV------------------------------------#

module load gcc/10.2.0 samtools/1.12 htslib/1.12
module load bwa/0.7.17

#--------------------INPUT----------------------------------#

echo "${SLURM_ARRAY_TASK_ID}"

INPUT=(`head -n ${SLURM_ARRAY_TASK_ID} ${1} | tail -n 1`) 

ID=$(  echo ${INPUT[@]} | cut -f 2 -d " ")  #e.g. 08UCPB10
R1=$(  echo ${INPUT[@]} | cut -f 4 -d " ")  #e.g. 08UCPB10_L001_R1.trimmed.fastq.gz
R2=$(  echo ${INPUT[@]} | cut -f 5 -d " ")  #e.g. 08UCPB10_L001_R2.trimmed.fastq.gz 
RGID=$(echo ${INPUT[@]} | cut -f 3 -d " ")  #e.g. 08UCPB10_L001
LB=$(  echo ${INPUT[@]} | cut -f 3 -d " ")  #e.g. 08UCPB10_L001
RG="@RG\tID:$RGID\tPL:ILLUMINA\tPU:unit1\tSM:$ID\tLB:$LB"

REF=/data/gpfs/projects/punim1823/genome/psco-genome.fasta.gz

# RG Definitions: 
# EXAMPLE: HER056678_S1_L001_R1_001.fastq
    # RGID* (read group ID; unique for each lane; this ID matches together PE reads: R1 and R2)=eg S1_L001
    # PL* (platform)=Illumina
    # PU=(use default)=unit1
    # SM* (sample name)=eg HER056678
    # LB* (DNA preparation library: for example this is "S")=eg S1


#--------------------OUTPUT--------------------------------#

outdir=/data/gpfs/projects/punim1823/snp-chip/agrf-illumina-data/fastp-results/bams

#--------------------RUN-----------------------------------#


# 1. Align with BWA and output as .bam (***note this method reqs more cores so may fail)
bwa mem -t 12 -R ${RG} ${REF} ${R1} ${R2} | samtools sort -o ${outdir}/${RGID}.sort.mapped.bam -

# 2. Index (creates CSI index)
samtools index -c ${outdir}/${RGID}.sort.mapped.bam