#!/bin/bash
#
#                             Remove duplicates
#  ***need to adjust mem and time reqs*****
#SBATCH -J Remove_Dups
#SBATCH --mem=50G  #Memory intensive: should be >40% more than required
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00:00  

#--------------------INFO-----------------------------------#

# Run from terminal as "sbatch merge-remDups.slurm merge-lists/<list.name> <individual.ID>""

#--------------------ERROR-CODES----------------------------#

set -o errexit
set -o pipefail
set -o nounset

#--------------------ENV------------------------------------#

module load gcc/10.2.0 samtools/1.12 htslib/1.12
module load picard/2.25.0-java-11

#--------------------INPUT----------------------------------#

LIST=$1  #list of multilibrary and/or multilane bam files with one file per line; e.g. Merge_lists/87930_mergelist.txt
ID=$2    #e.g. 87930

#--------------------PATHS----------------------------------#

picard=/usr/local/easybuild-2019/easybuild/software/core/picard/2.25.0-java-11/picard.jar
dir=/data/gpfs/projects/punim1823/snp-chip/agrf-illumina-data/fastp-results/bams
dir2=/data/gpfs/projects/punim1823/snp-chip/agrf-illumina-data/sandbox

#--------------------RUN-----------------------------------#

# 1. Merge multilibrary bam files 
samtools merge -b ${LIST} ${dir}/${ID}.merged.bam

# 2. Index 
samtools index -c ${dir}/${ID}.merged.bam

# 3. Remove duplicates with Picard tools 
java -jar -Xmx32G ${picard} MarkDuplicates \
        INPUT=${dir}/${ID}.merged.bam \
        OUTPUT=${dir}/${ID}.merged.nd.bam \
        METRICS_FILE=${dir2}/${ID}.metrics.txt TMP_DIR=/data/gpfs/projects/punim1823/snp-chip/agrf-illumina-data/tmp

# 4. Index 
samtools index -c ${dir}/${ID}.merged.nd.bam