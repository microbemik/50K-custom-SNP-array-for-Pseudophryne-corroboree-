# Mikaeylah Davidson, 2023 

# scripts used to design a 50K custom SNP array for Pseudophryne corroboree 

# all analysis was run on The Univeristy of Melbournes HPC 'SPARTAN'

# scripts are in numerical order 
# slurm jobs referenced throughout scripts are in alphabetical order 
    # these have been added for simplicity, but the alphanumeric need to be removed to run, or added to scripts

###--------------- SCRIPT 1 ---------------###
# Download psco genome and index 
1_Genome.sh
    # slurm scrips
A_samtools-index-file.slurm.sh
B_bwa-genome-index.slurm.sh


###--------------- SCRIPT 2 ---------------###
# Quality control and filtering of raw reads 
1_Quality-control.sh
    # slurm scrips
C_fastqc-run.slurm.sh
D_multiqc.slurm.sh
E_fastp-array.slurm.sh


###--------------- SCRIPT 3 ---------------###
# Align, merge and remove duplicates 
3_Align-merge.sh
    # slurm scrips
F_runBWA.slurm.sh
G_merge-remDups.slurm.sh


###--------------- SCRIPT 4 ---------------###
# Variant calling
4_Variant-calling.sh
    # slurm scrips
H_variantCalling-bcftools.slurm.sh


###--------------- SCRIPT 5 ---------------###
# Hard filtering
5_Hard-filtering.sh


###--------------- SCRIPT 6 ---------------###
# SNP filtering 
6_SNP-filtering.sh


###--------------- SCRIPT 7 ---------------###
# formatting for subimission
7_Submission-formating.sh
    # python script
I_71mer-format.py


