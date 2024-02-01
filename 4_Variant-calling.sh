# 
# Mikaeylah Davidson
# July 2023 
# PSCO custom SNP array 
# Variant calling
#

# using an adapted pipeline from
    # https://github.com/marine-omics/movp


###--------------- INPUT FILES ---------------###
# input file should be the merged, no duplicate, bam files
# eg. 
    B00699.merged.nd.bam



###--------------- CALL VARIANTS ---------------### 
# call SNPs using BCF tools 
    # bcftools/1.15

# run samples individually, and then merge
# create an array file with: 
    # line no, ID, file name, location
# list of IDs
awk -F "_" '{print $1}' ID | uniq > temp0
# array
awk '{print NR, $1, $1 ".merged.nd.bam " "/data/projects/punim1823/snp-chip/agrf-illumina-data/bams/" $1 ".merged.nd.bam" }' temp0 > bcftools-array

# run slurm 
sbatch variantCalling-bcftools.slurm bcftools-array

# output files = individual VCFs
# eg. 
    *.raw.bcftools.vcf.gz




###--------------- MERGE ---------------###
# individual read files are created, want one file with all reads
# merge file together using BCFtools 

# create a file list 
    # two lists, one for AGRF samples and one for VGP samples 
ls *bcftools* | awk '{print "/data/projects/punim1823/snp-chip/agrf-illumina-data/vcfs/" $1}' > list-agrf-bcftools
awk '{print "/data/projects/punim1823/snp-chip/vgp-illumina-data/vcfs/" $1 ".raw.bcftools.vcf.gz"}' IDs > list-bcftools-vgp

# files need to be indexed before merging: 
# load bcftools 
module load gcc/10.3.0 bcftools/1.15
# create list of files:
ls *bcftools* | awk '{print "bcftools index -c " $1}'
ls *bcftools* | awk '{print "bcftools index -c " $1}'

# copy and paste list into terminal, 
# eg.
#    bcftools index -c 08UCPB10.raw.bcftools.vcf.gz 

# merge AGRF & VGP samples
cat ~/projects/snp-chip/agrf-illumina-data/vcfs/list-bcftools-agrf ~/projects/snp-chip/vgp-illumina-data/vcfs/list-bcftools-vgp > list-bcftools-merged

# merge files 
    # use -l list-name
    #     -Oz = output compressed vcf
bcftools merge --file-list list-bcftools-merged -Oz -o psco.raw.bcftools.vcf.gz

# index merged file
bcftools index -c psco.raw.bcftools.vcf.gz 

# See summary stats of VCF file 
    # file is to large to view so saved to new file
bcftools stats psco.raw.bcftools.vcf.gz > bcftools-vcfOutput

# Summary = BCF Tools
    # number of samples	23
    # number of records	53,902,064
    # number of no-ALTs	0
    # number of SNPs	46,885,605
    # number of MNPs	0
    # number of indels	7,016,459
    # number of others	0
    # number of multiallelic sites	1,495,832
    # number of multiallelic SNP sites	183,620




###--------------- OUTPUT FILE  ---------------###
# output file is the merged SNP file, zipped VCF
# eg. 
psco.raw.bcftools.vcf.gz
