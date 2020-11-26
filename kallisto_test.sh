#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


# 23/10/2020: KALLISTO

RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
KALLISTO=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/Kallisto
All_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
TG_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/TG_Merged
WT_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/WT_Merged
FLAIR_OUTPUT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/FLAIR

module load Miniconda2
run_kallisto(){
    source activate sqanti2
    
    echo "Processing Kallisto for $1"
    cd $4
    kallisto version
    time kallisto index -i $1_Kallisto.idx $2 2> $1_Kallisto.index.log
    time kallisto quant -i $4/$1_Kallisto.idx --fr-stranded $3/$1_R1.fq --rf-stranded $3/$1_R2.fq -o $4 2> $1_Kallisto.quant.log
    mv abundance.tsv $1".abundance.tsv"
    
    # problem: retained co-ordinates, which does not input well into SQANTI2
    echo "Kallisto original $1.abundance.tsv"
    head $1".abundance.tsv"
    # solution: retain the PB.ID
    while read line ; do
      first=$( echo "$line" |cut -d\| -f1 ) # each read line, and remove the first part i.e. PacBio ID
      rest=$( echo "$line" | cut -d$'\t' -f2-5 ) #save the remaining columns
      echo $first $rest # concatenate
    done < $4/$1".abundance.tsv" > $4/$1.temp.abundance.tsv
    
    header=$(head -n 1 $4/$1.abundance.tsv)
    sed -i '1d' $4/$1.temp.abundance.tsv # remove header of temp.file to be replaced
    echo $header > foo
    cat foo $4/$1.temp.abundance.tsv > $4/$1.mod.abundance.tsv 
    
    echo "Kallisto $1.mod.abundance.tsv"
    head $4/$1.mod.abundance.tsv
    rm $1.temp.abundance.tsv
    rm foo

    source deactivate
}

#run_kallisto All_RNASeq $All_Isoseq3_WKD/TOFU/mm10/All_Merged.collapsed.filtered.rep.fa $KALLISTO $KALLISTO
#run_kallisto TG_RNASeq $TG_Isoseq3_WKD/TOFU/mm10/TG_Merged.collapsed.filtered.rep.fa $KALLISTO $KALLISTO
#run_kallisto WT_RNASeq $WT_Isoseq3_WKD/TOFU/mm10/WT_Merged.collapsed.filtered.rep.fa $KALLISTO $KALLISTO
run_kallisto All_RNASeq $FLAIR_OUTPUT/flair.collapse.isoforms.fa $KALLISTO $FLAIR_OUTPUT
