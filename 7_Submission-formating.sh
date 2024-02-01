# 
# Mikaeylah Davidson
# Aug 2023 
# PSCO custom SNP array 
# Formatting for submission
#


###--------------- INPUT FILES ---------------###
# input file is the completely filtered SNP file
SNPs-thin1kb.recode.vcf
#  no. SNPs = 537,368 


###--------------- SUBMISSION PARAMETERS ---------------###
# 
# file MUST be submitted to specific requirments set by ThermoFisher 
# need file in 71-mer format 
# 

# extract snp with flanking sequence - 71mer format 
# 35bp[snpRef/snpAlt]35bp 
# following: https://www.biostars.org/p/334253/
    # finswimmers python script 
    # script eddited to add species, ID, importance and the format.
        # script saved as == 71mer-format.py
# 
# run python file 
python 71mer-format.py > psco_SNPs_537368_71merFormat.txt

# wc
537,368

# add header to file 
sed -i $'1 i\\\torganism\tsnpid\tseventyonemer\timportance\tchromosome\tposition\tref_allele\talt_allele' psco_SNPs_537368_71merFormat.txt

# check columns are correct, and values line up
awk '{print $1}' psco_SNPs_537368_71merFormat.txt | head 

 