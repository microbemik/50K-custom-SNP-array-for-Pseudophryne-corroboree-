#
# Python script to take VCF file of SNPs, and convert to 71mer format
# needed for SNP submission to Thermofisher 
# 
# From: https://www.biostars.org/p/334253/
# following finswimmers post
# 
# script edited to add in columns required for submittion, and formatted
# 
# load modules first:
# module load GCC/11.3.0 Pysam/0.19.1
# 

import pysam


# open vcf file
vcf = pysam.VariantFile("snpFilt-step4.1.vcf")
# open fasta file
genome = pysam.FastaFile("psco-genome-cp.fasta")
# define by how many bases the variant should be flanked
flank = 35

# iterate over each variant
for record in vcf:
    # extract sequence
    #
    # The start position is calculated by subtract the number of bases
    # given by 'flank' from the variant position. The position in the vcf file
    # is 1-based. pysam's fetch() expected 0-base coordinate. That's why we
    # need to subtract on more base.
    #
    # The end position is calculated by adding the number of bases
    # given by 'flank' to the variant position. We also need to add the length
    # of the REF value and subtract again 1 due to the 0-based/1-based thing.
    #
    # Now we have the complete sequence like this:
    # [number of bases given by flank]+REF+[number of bases given by flank]
    seq = genome.fetch(record.chrom, record.pos-1-flank, record.pos-1+len(record.ref)+flank)

    # print out tab seperated columns:
    # CRHOM, POS, REF, ALT, flanking sequencing with variant given in the format '[REF/ALT]'
    print(
        "Pseudophryne_corroboree",
        "psco_" + record.id,
        '{}[{}/{}]{}'.format(seq[:flank], record.ref, record.alts[0], seq[flank+len(record.ref):]),
        "Standard",
        record.chrom,
        record.pos,
        record.ref,
        record.alts[0],
        sep="\t"
    )