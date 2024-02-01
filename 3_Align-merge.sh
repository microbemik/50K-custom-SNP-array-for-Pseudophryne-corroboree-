# 
# Mikaeylah Davidson
# June 2023 
# PSCO custom SNP array 
# Align & merge files 
#




###--------------- INPUT FILES ---------------###
# input files are the trimmed files, they should be zipped fastq 
# eg. 
    B00699_L003_R1.trimmed.fastq.gz
#   ID, Lane, F/R





###--------------- BWA ---------------### 
# align read files to the genome using BWA 

# create array for BWA
# print location of trimmed files
ls `pwd`/*trimmed.fastq.gz | sort > trimmed-fastq.list

# want ID, and then unique ID and lane ID 
    # eg. 9UO05 9UO05_L001
ls *gz | awk -F "_" '{print $1, $1"_"$2}' > temp-1
    # this list doubled up as both R1 and R2 so total = 156 IDs 

# unique IDs only 
sort temp-1 | uniq | wc 
    # wc = 78 unique file names
# over ride IDs file 
# print row number 
sort temp-1 | uniq | awk '{print NR, $1, $2}' > IDs

# remove temp-1 file
rm temp-1

# pull R1 and R2 and pathw
grep "R1" trimmed-fastq.list > R1
grep "R2" trimmed-fastq.list > R2

# combine  
paste -d " " IDs R1 R2 
# make into array file 
    # == bwa-array
paste -d " " IDs R1 R2 | awk '{print $0}' > bwa-array

# run BWA slurm 
sbatch runBWA.slurm bwa-array
    # output files = *sort.mapped.bam




###--------------- MERGE & REMOVE DUPLICATES ---------------###

# Merge list
# create a list of BAM the files which need to be merged
# each animal has 2-4 sample runs, so 2-4 BAM files
# create a list for each animal of the files which need to be merged 
ls *R1*gz | wc
    # 78 files 
# want file name and lane with .sort.mapped.bam
ls *R1*gz | awk -F "_" '{print $1"_"$2".sort.mapped.bam"}' > temp1
# add path to file 
ls *R1*gz | awk -F "_" '{print "/data/gpfs/projects/punim1823/snp-chip/agrf-illumina-data/fastp-results/bams/" $1"_"$2".sort.mapped.bam"}' > temp1
    # eg. output
        # /data/gpfs/projects/punim1823/snp-chip/agrf-illumina-data/fastp-results/bams/B40307_L004.sort.mapped.bam

# want a grep command, which will use the unique ID to pull each lane find into a merge list
# create a temp file with the command, which can then copy and pasted into terminal
# create list of unique IDs
    # already have a file of IDs called "IDs"
awk '{print $2}' IDs | sort | uniq > temp2    

# create file of grep command 
awk '{print "grep " "\"" $1 "\"" " temp1" " > " $1 "_mergelist.txt"}' temp2 > temp3

# pull list of files into the unique merge list
# open temp3
# copy and paste contents of file into terminal and run
# should make *_mergelist.txt files, with the different bam files for each ID in it
# check files are correct
    # eg. 
    #  cat B00699_mergelist.txt

# run slurm script to merge 
sbatch merge-remDups.slurm merge-lists/<list.name> <individual.ID>




###--------------- OUTPUT FILE  ---------------###
# output file should be the merged, no duplicate, bam files
# eg. 
B00699.merged.nd.bam

