# 
# Mikaeylah Davidson
# July 2023 
# PSCO custom SNP array 
# Hard filtering 
#

# using an adapted pipeline from
    # https://github.com/marine-omics/movp
    # https://speciationgenomics.github.io/filtering_vcfs/


###--------------- INPUT FILES ---------------###
# input file is the merged SNP file, zipped VCF
psco.raw.bcftools.vcf.gz




###--------------- REQUIREMENTS ---------------###
# Hard filtering requires two things: 
#   1. List of bcftools-generated vcf files to merge
        # done in 4_Variant-calling.sh
#   2. BED file with list of CHR BEG END to keep. 
        # create from genome index
grep -v "JAQ" psco-genome.fasta.fai | awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' > keep-regions.bed




###--------------- FILTERING PARAMETERS ---------------###
# STEPS
1. Remove unasigned contigs
2. Sequence quality = QUAL >20 
3. Read depth = MIN(DP)>5 & MAX(DP)<25
4. Remove individuals with >70% missingness 
5. Missing data per SNP 20% = F_MISSING<0.2

# set variables
FILE=vcfs/psco.raw.bcftools.vcf.gz
LIST2=/data/projects/punim1823/snp-chip/merged/keep-regions.bed  #bed file with list of CHR BEG END to keep
dir=/data/projects/punim1823/snp-chip/merged/vcfs


###--------------- No. SNPs ---------------###
# check number of SNPs begining with in raw unfiltered file
bcftools stats psco.raw.bcftools.vcf.gz | head -n 50
# 46,885,605 SNPs


###--------------- HARD FILTERING ---------------###
# STEP 1
# Pull SNPs and keep only chrs
bcftools view ${dir}/psco.raw.bcftools.vcf.gz -v snps --regions-file ${LIST2} -Oz -o ${dir}/file1.bcftools.vcf.gz
# Update header (removed uncharaterized contigs)
# extract header
bcftools head ${dir}/file1.bcftools.vcf.gz > ${dir}/vcf-header
# remove JAQ rows
grep -v "JAQ" ${dir}/vcf-header > ${dir}/vcf-header2
# replace header
bcftools reheader -h ${dir}/vcf-header2 ${dir}/file1.bcftools.vcf.gz > ${dir}/file2.bcftools.vcf.gz
# index
bcftools index -c ${dir}/file2.bcftools.vcf.gz

# results 
bcftools stats file2.bcftools.vcf.gz | head -n 50
        #  no. SNPs              = 45,579,322
        #  no. SNPs removed      = 1,306,283 



# STEP 2
# sequences QUAL > 20
bcftools view ${dir}/file2.bcftools.vcf.gz -i 'QUAL>20' -Oz -o ${dir}/file3.bcftools.vcf.gz

# results 
bcftools stats file3.bcftools.vcf.gz | head -n 50
        #  no. SNPs              = 37,379,322
        #  no. SNPs removed      = 8,200,000 


# STEP 3
# read depth = MIN(DP)>5 & MAX(DP)<25
        # filter INFO/DP : read depth per sample. So if what to cut off 5<DP<25 and had 52 samples. WE should multipy number of DP with number of sample
        #   For example: e 'FMT/DP<5 | FMT/DP>25'  -->  e 'FMT/DP<250 | FMT/DP>1250' 
# minimum depth
bcftools view ${dir}/file3.bcftools.vcf.gz -i 'INFO/DP>115' -Oz -o ${dir}/file3.1.bcftools.vcf.gz
# maximum depth
bcftools view ${dir}/file3.1.bcftools.vcf.gz -i 'INFO/DP<575' -Oz -o ${dir}/file4.bcftools.vcf.gz

# results 
bcftools stats file4.bcftools.vcf.gz | head -n 50
        #  no. SNPs              = 15,763,936
        #  no. SNPs removed      = 21,615,386


# STEP 4
# remove individuals with >70% missingness 
# generate missingness per inviduals 
vcftools --gzvcf file4.bcftools.vcf.gz --missing-indv --out missingness_report
# F_MISS = the fraction of missing genotypes for each individual 
            # It is calculated as the ratio of the number of missing genotypes to the total number of genotypes for that individual.
        # eg. if an individual has F_MISS value of 0.4, it means that 40% of their genotypes are missing in the VCF file.
# NO indvidual has >70% missingness 
# range is 32-44%



# STEP 5
# Missing data per SNP 20% = F_MISSING<0.2
    # 3 animals can have missing genotypes per SNP
bcftools view file4.bcftools.vcf.gz -i 'F_MISSING<0.2' -Oz -o file5.bcftools.vcf.gz

# results 
bcftools stats file5.bcftools.vcf.gz | head -n 50
        #  no. SNPs              = 3,280,265
        #  no. SNPs removed      = 12,483,671 




###--------------- UPDATE HEADER ---------------###
# need to annotate labels = 
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' file5.bcftools.vcf.gz -O z -o bcftools.SNPs.vcf.gz

# check
bcftools view -H bcftools.SNPs.vcf.gz | wc -l
# 3,280,265 SNPs

#  check the number of variants per chromosome
bcftools query -f '%CHROM\t%POS\n' bcftools.SNPs.vcf.gz > chr-pos
awk '{print $1}' chr-pos | sort | uniq -c | awk '{print $2, $1}' | sort -V > variants-per-chr
cat variants-per-chr
        # 1 534798
        # 2 411944
        # 3 304133
        # 4 344549
        # 5 380754
        # 6 340779
        # 7 240331
        # 8 189547
        # 9 166408
        # 10 162744
        # 11 117539
        # 12 86739




###--------------- OUTPUT FILE  ---------------###
# output file is the hard filtered SNP file, zipped VCF
bcftools.SNPs.vcf.gz
#  no. SNPs = 3,280,265