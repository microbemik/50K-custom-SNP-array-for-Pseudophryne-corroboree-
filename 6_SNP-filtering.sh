# 
# Mikaeylah Davidson
# Aug 2023 
# PSCO custom SNP array 
# SNP filtering 
#


###--------------- INPUT FILES ---------------###
# input file is the hard filtered SNP file, zipped VCF
bcftools.SNPs.vcf.gz
#  no. SNPs = 3,280,265




###--------------- FILTERING PARAMETERS ---------------###
# STEPS
1. Remove low-complexity regions
2. Remove multiallelic, monomorphic, and SNPs within 100bp of indels
3. Filter HWE variance with p<10^6
4. Filter MAF 0.1
5. Check for duplicated SNPs
6. Remove SNP wildcards - “N” errors
7. Remove SNP which have any unfiltered SNP within 31bp either side
8. Thining by distance 

###--------------- No. SNPs ---------------###
# check number of SNPs begining with in raw unfiltered file
bcftools bcftools.SNPs.vcf.gz | head -n 50
# 3,280,265 SNPs


###--------------- STEP 1: LOW COMPLEXITY REGIONS ---------------###
# 
# flag regions of low complexity in genome and remove 
# 


# use DUSTMASKER
    # https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/README
# 
# install using conda 
# load anaconda
module load anaconda3/2021.11
#making env instead
conda create -name blast 
# activate environment 
eval "$(conda shell.bash hook)"
conda activate blast
    # deactivate enviro
    #    conda deactivate
# instal 
conda install -c bioconda blast
#Check installation
dustmasker -h

# create list of low complexity regions in genome
# 
# use unzipped, genome in fasta format
dustmasker -in psco-genome.fasta -infmt fasta -outfmt acclist -out psco-genome-dustmasker.acclist
    # eg of list:
        #chromosome    #start.pos  #end.pos
            # >1      1641          1697
            # >1      1765          1771
            # >1      2061          2067
# count number of positions = no.lines
wc -l psco-genome-dustmasker.acclist
    # == 12,806,735       #includes header 

# change list into bed format 
cp psco-genome-dustmasker.acclist psco-genome-dustmasker-acclist.bed

# remove the ">" from infront of the chromosome number 
# 
# ***IF THE CHROMOSOME NAMES ARE NOT EXACTLY THE SAME IT WILL NOT WORK!!!!!*** 
# 
# remove ">"
sed 's/^.//' psco-genome-dustmasker-acclist.bed > psco-genome-LClist.bed
# check 
head psco-genome-LClist.bed 
    # ">" removed

# generate list of SNPs within the low complexity regions
# 
# Find common SNP between LC list file and vcf files 
# using bedtools intersect
    # VCF file = bcftools.SNPs.vcf.gz
bedtools intersect -b ../hardFilt/bcftools.SNPs.vcf.gz -a ../psco-genome-LClist.bed -wb > lowcomplex-regions-bedtools.bed     
    # give all SNPs in VCF which overlap with lowcomplexity regions in LC list file
# 91,484 SNPs in LC regions  
    # file == lowcomplex-regions-bedtools.bed 
# use to remove snps which fall into these regions from the SNP file 


# remove SNPs within the low complexity regions
# 
# check no of SNPs 
bcftools view -H ../hardFilt/bcftools.SNPs.vcf.gz | wc -l  
    # == 3,280,265 SNPs

# remove SNPs in low complexity regions 
bcftools view -T ^lowcomplex-regions-bedtools.bed ../hardFilt/bcftools.SNPs.vcf.gz -Oz > snpFilt-step1.vcf.gz 

# results
bcftools view -H snpFilt-step1.vcf.gz  | wc -l     
    # no SNPs in new file == 3,188,781
    # no SNPs in orginal file == 3,280,265
    # no SNPs removed == 91,484
# == matches the number of snps in the LC file 
# completed. 




###--------------- STEP 2: multiallelic, monomorphic and indel regions ---------------###
# 
# remove multiallelic SNPs, monomorphic SNPs, and SNPs within 100bp of indels
# 


# exclude -e all sites at which no alternative alleles are called for any of the samples 
#   AC==0 
# and all sites at which only alternative alleles are called 
#   AC==AN
# remove SNPs that are within 100 bp of indels, we add the flag 
    # --SnpGap 100
# want only biallelic SNPs, we use the 
#   bcftools view -m2 -M2 -v snps
# which set both the minimum -m and maximum -M number of alleles to 2, and keeps only snps via -v snps
# use -e and the logical operator “or” || to exclude low-quality variants from the file
bcftools filter -e 'AC==0 || AC==AN' --SnpGap 100 snpFilt-step1.vcf.gz | bcftools view -m2 -M2 -v snps -O z -o snpFilt-step2.vcf.gz

# results
bcftools view -H snpFilt-step2.vcf.gz | wc -l   
    # no SNPs == 3,146,940
    # no SNPs removed == 41,841




###--------------- STEP 3: FILTER HWE ---------------###
# 
# filter Hardy-Weinberg equilibrium 10^6 
# 

vcftools --gzvcf snpFilt-step2.vcf.gz --hwe 0.000001 --recode --recode-INFO-all --out snpFilt-step3

# results 
bcftools view -H snpFilt-step3.recode.vcf | wc -l  
    # no of SNPs == 3,146,940
    # no SNPs removed == 0

# this popuation should not be in HWE
# its predominatly males, from geographical separate populations 
# they were the last wild caught animals from this species 




###--------------- STEP 4: FILTER MAF ---------------###
# 
# filter minor allel frequency 0.1
#

vcftools --vcf snpFilt-step3.recode.vcf  --maf 0.1 --recode --recode-INFO-all --out snpFilt-step3.1  

# results 
bcftools view -H snpFilt-step3.1.recode.vcf | wc -l  
    # no of SNPs == 2,798,218
    # no SNPs removed == 348,722




###--------------- STEP 5: REMOVE DUPLICATES ---------------###
# 
# remove SNPs which are duplicated 
#


# convert vcf to bed & mark snps 
plink --vcf snpFilt-step3.1.recode.vcf --make-bed --out snpFilt-step4     
# convert plink.bed file to regular .bed file 
plink --bfile snpFilt-step4 --recode vcf bgz --out snpFilt-step4.0
# use bedops to conver to UCSC .bed
    # installed bedops in conda environment: 
zcat snpFilt-step4.0.vcf.gz | vcf2bed - > snpFilt-step4.0.bed

# Extrac SNP position using awk
awk '{print $1,$2,$3;}' OFS='\t' snpFilt-step4.0.bed > snp-positions.txt
# This step takes 150bp region eith side of the SNP position 
# Note when using awk the space may not be tab delimited, if its not the file can't be read
    # use OFS='\t' to ensure it is
awk '{$2-=150;$3+=150};{print}' OFS='\t' snp-positions.txt > snp-positions-150bp.bed

# check no. SNPs is still matching 
wc -l snp-positions-150bp.bed 
# 2,798,218   == yes

# convert .bed to .fasta
# this keeps the bed regions, and pulls the sequences from those regionins in the genome into the file 
bedtools getfasta -fi psco-genome-cp.fasta -bed snp-positions-150bp.bed -fo snp-positions-150bp.fasta
# the number of lines in this file should double, compared to the .bed file, as there is now 1 line of sequence 


# map the sequences to the genome 
# using BWA
#  The BWA-MEM algorithm performs local alignment
## Code: $ bwa mem index_prefix [input_reads.fastq|input_reads_pair_1.fastq input_reads_pair_2.fastq] [options]
#  index_prefix is the index for the reference genome generated from bwa index,
#  input_reads.fastq, input_reads_pair_1.fastq, input_reads_pair_2.fastq are the input files of sequencing data that can be single-end or paired-end respectively
#  Additional options for bwa mem can be found in the BWA manual.
#  bwa mem -t 10 index_prefix input_reads_pair_1.fastq input_reads_pair_2.fastq > bwa_mem_alignments.bwa.sam
# bwa mem ref.fa reads.fq > aln-se.sam
bwa mem -t 10 psco-genome-cp.fasta snp-positions-150bp.fasta > snp-150bp-mappedgenome.sam
# this likely will need to be submitted as slurm job as quite intensive

# convert sam to bam
samtools view -b snp-150bp-mappedgenome.sam > snp-150bp-mappedgenome.bam

# print the first column to find dupilcated values 
#  
samtools view snp-150bp-mappedgenome.sam | awk '{print $1}' | head
samtools view snp-150bp-mappedgenome.sam | awk '{print $1}' > snp-150bp-mapped-dups.txt

# Find duplicated   
awk 'seen[$0]++ == 1' snp-150bp-mapped-dups.txt 
# no results

wc -l snp-150bp-mapped-dups.txt 
# no SNPs 2,798,218
# no duplicated sequences 


###--------------- STEP 6: REMOVE WILDCARDS ---------------###
# 
# identify sequences contain sequence wildcard "N" and remove  
#

awk '$10 ~ "N" {print $1,$10}' OFS='\t' snp-150bp-mappedgenome.sam > snp-150bp-mappedgenome-N.txt 

wc -l snp-150bp-mappedgenome-N.txt
# 92 SNPs contain N 


# create a file of the snp positions, from the regions 
# use grep to remove
awk '{print $1}' snp-150bp-mappedgenome-N.txt > snp-150bp-toremove.txt
awk -F "[:|-]" '{print $2}' snp-150bp-toremove.txt > snp-position150bp-toremove.txt
# as these values are 150 bp from the target snp, need to get the actual snp position 
    # add 151 to the snp poistion 
awk '{$1+=151};{print}' snp-position150bp-toremove.txt > snp-toremove.txt
# 92 snp positions to remove: 
    
# create grep command, to be able to copy and paste the 92 snps to remove 
# pipe all the snp positions together
    # use -w flag so it only pulls the perfect match
# unzip file 
gzip -d snpFilt-step4.0.vcf.gz
grep -vwf snp-toremove.txt snpFilt-step4.0.vcf > snpFilt-step4.1.vcf

# check that only 92  snps have been removed 
bcftools view -H snpFilt-step4.1.vcf | wc -l  
# no SNPs == 2,798,126   = yes only 92 removed




###--------------- STEP 7: FLANKING REGIONS ---------------###
# 
# Remove SNP which have any unfiltered SNP within the flanking region = 31bp either side
# 31 is used as it counts the SNP as pos#1 
# 

# need:
# 1. file of unfiltered SNPs == file3.bcftools.vcf.gz
    # using file which has had SNPs in unassigned contigs removed, and quality<20 removed
# 2. current filtered snp file == snpFilt-step4.1.vcf


# check no. reccords/SNPs
bcftools view -H file3.bcftools.vcf.gz | wc -l 
# 37,379,322 snps
bcftools view -H snpFilt-step4.1.vcf | wc -l 
# 2,798,126 snps 

# 37,379,322 - 2,798,126 = 34,581,196 
#  === want a file with 34,581,196 records 


# Step 1: 
# remove filtered snps fropm orginal snp file 
# zip filtered file 
bgzip -c snpFilt-step4.1.vcf > snpFilt-step4.1.vcf.gz

# using bcftools: 
    # https://www.biostars.org/p/392324/
bcftools view -T ^snpFilt-step4.1.vcf.gz file3.bcftools.vcf.gz > filtSNPs-subtracted.vcf
bcftools view -H filtSNPs-subtracted.vcf| wc -l 
# 34,581,196
# removed correct amount 

# Step 2: 
# create bed file with flanking regions around snps on original file 
# 31bp flanking region 

# first convert vcf to bed, 
# zip file 
bgzip -c filtSNPs-subtracted.vcf > filtSNPs-subtracted.vcf.gz
# and then use the bed file to add flanking region 
zcat filtSNPs-subtracted.vcf.gz | vcf2bed - > filtSNPs-subtracted.bed

# then add 31bp flanking region to either side of SNP
# NOTE - 
    # this splits the forward and reverse flanking region across two lines 
    # first the 5" flank, then the 3" flank
bedtools flank -i filtSNPs-subtracted.bed -g psco-genome-cp.fasta.fai -b 31 > 31bp-flanking.bed


# Step 3: 
# intersect filtered snp file with orig snp bed file & 
# remove any filtered snps which fall in the regions surrounding orignal file 

# use bedtools intersect 
    # https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
# -v reports the absence of any overlapping features
# which “A” features do not overlap with any of the “B” features
# eg. which SNPs don’t overlap with 30bp either side of SNPs from original file 
bedtools intersect -a snpFilt-step4.1.vcf.gz -b 31bp-flanking.bed -v > flankingSNPs-removed.vcf
# run as a slurm script

wc -l flankingSNPs-removed.vcf
# no SNPs == 2,222,878  
# no SNPs removed == 575,248




###--------------- STEP 8: THINNING ---------------###
# 
# thinning SNPs by distance 
# aiming to submit ~ 500,000 SNPs for evaluation to Thermofisher   
#

# Thinning 1kb
vcftools --vcf ../flankingSNPs-removed.vcf --recode --recode-INFO-all --thin 1000 --out SNPs-thin1kb

# results 
bcftools view -H SNPs-thin1kb.recode.vcf   | wc -l 
# no SNPs == 537,368 




###--------------- OUTPUT FILE  ---------------###
# output file is the filtered SNP file, not zipped VCF 
SNPs-thin1kb.recode.vcf
#  no. SNPs = 537,368 

