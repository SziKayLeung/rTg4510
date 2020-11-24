#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrcq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

module load R/3.4.3-foss-2016b-X11-20160819
module load Miniconda2
source activate nanopore
# 24/07/2019: conda install -c bioconda blast

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
AARON_REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/reference
TARGETED_HUMAN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/human

### Prepared Probes for BLAST
# Loaded SziKay_Human23072019.txt to $TARGETED_HUMAN
echo "Processing Human Probes (23/07/2019)  for blast"
# Process txt file to include ">" at the beginning of each row, column 1, extract only first two columns, remove first two header rows
sed 's/^/>/' SziKay_Human_23072019.txt | awk '{print $1"\n"$2}' | sed '1,2d' > blast_Human_input.txt
head blast_Human_input.txt

### Prepared Reference Genome for BLAST
# Copied Aaron's hg38 to from his reference folder to $REFERENCE
# cp $AARON_REFERENCE/hg38.fa $REFERENCE
cd $TARGETED_HUMAN
makeblastdb -in $REFERENCE/hg38.fa -dbtype nucl -out hg38.blast

### Blast Probes to Reference Genome
echo "Blast Human Probes to hg38"
blastn -query $TARGETED_HUMAN/blast_Human_input.txt -db hg38.blast -out blast_Human_output.txt -outfmt 6
head blast_Human_output.txt

### PreProcess Blasted Probes 
# some probes have START and END coordinates switched around, therefore not recognised by USCS Table Browser
# Rscript Switch$Prepare_Probes.r <path/blast_output_file> <path/output_directory>
Rscript Switch&Prepare_Probes.r $TARGETED_HUMAN/blast_Human_output.txt $TARGETED_HUMAN

### Convert Blast output coordinates to Gene Names on USCS Table Browser 
# download blast_Human_output.txt 
# on USCS table browser: https://genome.ucsc.edu/cgi-bin/hgTables
    # clade: Mammal; genome: Human; assembly: hg38 
    # track: NCBI RefSeq 
    # click defined regions, change, and upload blast_Human_output.txt 
    # output_format: selected fields from primary and related tables 
    # get output
    # select "name", "chrom", "score"
    # save output as "Blast_Hg38_output.txt"

### Find number of hits from Blast_Hg_38_output.txt
cut -f4 Blast_Hg38_output.txt | sed '1d' | sort | uniq > genes_blast.output.txt

while read line
do
   counts=$(grep $line Blast_Hg38_output.txt | wc -l)
   echo $line
   echo $counts
done < genes_blast.output.txt >> gene_blast_hits.txt

while read line
do
   counts=$(grep $line Blast_Hg38_output.txt | wc -l)
   echo $line
   #echo $counts
done < genes_blast.output.txt