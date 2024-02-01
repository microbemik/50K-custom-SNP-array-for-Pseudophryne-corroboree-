# 
# Mikaeylah Davidson
# April 2023 
# PSCO custom SNP array 
# Genome
#




###--------------- GENOME ---------------### 
# in order to align files, need genome and index files

# download psco genome from ncbi
    https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Pseudophryne_corroboree/latest_assembly_versions/GCA_028390025.1_aPseCor3.hap2/
        # genome file =
            GCA_028390025.1_aPseCor3.hap2_genomic.fna.gz    
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Pseudophryne_corroboree/latest_assembly_versions/GCA_028390025.1_aPseCor3.hap2/GCA_028390025.1_aPseCor3.hap2_genomic.fna.gz
# rename file 
mv GCA_028390025.1_aPseCor3.hap2_genomic.fna.gz genome-original.fna.gz
    # new name = genome-original.fna.gz

# count number of base pairs 
zgrep -v ">" genome-original.fna.gz | wc
    # === 110910913 110910913 8983669706
        # lines     words     characters
    # to get the no. of bp in the file == characters - lines
zgrep -v ">" genome-original.fna.gz | wc | awk '{print $3-$1}'
    # === 8872758793 bp 




###--------------- RENAME CHROMOSOMES ---------------###

zcat genome-original.fna.gz | head
    # names all start with ">"

# view first 20 chromosomes 
zgrep ">" genome-original.fna.gz | head -n20
    # chromosomes apear to start with ">CM0"
        # ">JAQOLT" are unasigned contigs within the chr

# veiw only chromosomes = should be 12
zgrep ">CM0" genome-original.fna.gz
    # yes = 12

# pull chromosome names into new file 
zgrep ">CM0" genome-original.fna.gz > ncbi-chrs-names
    # file == "ncbi-chrs-names"

# check file length
wc -l ncbi-chrs-names
    #  === 12

# check file order 
cat ncbi-chrs-names
    # === files in numerical order 1-12

# make name replacement file 
    # remove ">" = replace with ' '
    # \t = replaces spaces with tabs
sed 's/>//g' ncbi-chrs-names | awk 'OFS="\t" {print $0, NR}' > new-chrs-names
    # new file == 'new-chrs-names'

# create copy of genome 
cp genome-original.fna.gz genome-copy.fna.gz

# rename chromosomes in copy file 
    # want chromosomes to be named "1", "2".... etc
zcat genome-copy.fna.gz | awk -F "\t" 'NR==FNR{array[">"$1]=(">" $2);next} {if($1 in array){$1=array[$1]}}1' new-chrs-names - > genome-chr-renamed.fasta
    # new file == "genome-chr-renamed.fasta"

# check format of new.fasta file 
head genome-chr-renamed.fasta

# bgzip file
bgzip genome-chr-renamed.fasta

# rename file 
mv genome-chr-renamed.fasta.gz psco-genome.fasta.gz
    # genome with renamed chromosomes === "psco-genome.fasta.gz"




###--------------- INDEX GENOME ---------------###

# create a FAIDX index
    # using samtools/1.12 

# run using slurm script
sbatch samtools-index-file.slurm

#  check file 
cat psco-genome.fasta.gz.fai | head -n20

#  pull sizes
cut -f1,2 psco-genome.fasta.gz.fai > chrs-sizes

# only want chromosomes, not unassigned contigs
grep -v "J" chrs-sizes > chrs-sizes-assigned

# summarise the % of genome and chromosome
FILE=psco-genome.fasta.gz.fai
LEN=$(cut -f1,2 ${FILE} |  awk '{sum+=$2} END {print sum}')
CHR=$(grep -v "JAQ" ${FILE} | cut -f1,2 | awk '{sum+=$2} END {print sum}')
awk -v LEN="$LEN" -v CHR="$CHR" 'BEGIN {print CHR/LEN*100}'
    # === 92.3081 




# index whole genome
    # using BWA 

# run using slurm script
sbatch bwa-genome-index.slurm



# NOTE:
# may need an unzipped genome indexed for some purposes
# use same scripts, but unzip genome and run. 
