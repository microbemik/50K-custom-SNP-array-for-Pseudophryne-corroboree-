# 
# Mikaeylah Davidson
# April 2023 
# PSCO custom SNP array 
# Quality control & filtering 
#


###--------------- INPUT FILES ---------------###
# input files are the raw read files, they should be zipped fastq 
# eg. 
    A80101_L001_R1.fastq.gz
#   ID, Lane, F/R




###--------------- FASTQC ---------------### 
# check quality of data using FastQC 
    # http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# create array list of file names, with line numbers into new file
ls *gz | awk '{print NR, $0}' > array-fastqc

# load module 
# fastqc/0.11.9-java-11
module spider fastqc 

# run slurm script
sbatch fastqc-run.slurm fastqc-array

# run MultiQC to combine all FastQC 
    # https://multiqc.info/docs/getting_started/running_multiqc/

# load module 
module spider multiqc

# run slurm script 
sbatch multiqc.slurm




###--------------- FASTP ---------------###
# trimm sequences using FASTP
            #  https://speciationgenomics.github.io/fastp/#:~:text=These%20polyG%20tails%20need%20to,quality%20parts%20are%20trimmed%20off.
# removes
    # adaptors 
    # polyG tails 
        # due to using NOVASEQ sequencing, need to remove polyG tails 
    # -l 50 = removed sequences under 50bp after trimming 
    # uncalled or unqualified bases 

# load module 
# fastp/0.23.2
module spider fastp

# create array with R1 and R2 files names on the same line 
ls *R1*gz | awk '{print $0}' > R1
ls *R2*gz | awk '{print $0}' > R2
# also want the animal ID seperate to be able to re-name the output files 
    # want all information before the R1/2
awk -F "_" '{print $1"_"$2}' R1 > ID

# merge files together into the array 
paste -d " " ID R1 R2 | awk '{print NR, $0}' > fastp-array

# run slurm script 
sbatch fastp-array.slurm fastp-array


###--------------- OUTPUT FILE  ---------------###
# output file should be the trimmed raw reads file, they should be zipped fastq  
# eg. 
B00699_L003_R1.trimmed.fastq.gz

