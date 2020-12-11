#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --error=lengths.error
#SBATCH --output=lengths.out

FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq
CCS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/Isoseq3_WKD/CCS
LENGTHS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/LENGTHS
SAMPLES=(Q21 K17 M21 O23 S23 K23 O18 K18 S18 L22 Q20 K24)
SAMPLES_MERGED=(All_Merged TG_Merged WT_Merged)
RAW_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/rawdata


#clustered  
INDIVIDUAL_CLUSTERED=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/Isoseq3_WKD/CLUSTER

TG4510=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510
MERGED_CLUSTERED=($TG4510/All_Merged/CLUSTER $TG4510/TG_Merged/CLUSTER $TG4510/WT_Merged/CLUSTER)

# tofu
INDIVIDUAL_TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/TOFU
MERGED_TOFU=($TG4510/All_Merged/TOFU/mm10 $TG4510/TG_Merged/TOFU/mm10 $TG4510/WT_Merged/TOFU/mm10)

source $FUNCTIONS/Post_Isoseq3_Functions.sh

# convert_bam2fq <input/output_prefix_name> <input_dir> <output_dir>
# extract_length_fastq <input_file> <output_prefix> <output_dir>
#for i in ${SAMPLES[@]}; do convert_bam2fq $i "_ccs" $CCS $CCS; done
# convert_bam2fq K17 $CCS $LENGTHS
#extract_length_fastq $SAMPLES.fastq $SAMPLES $CCS 

# individual samples
#for i in ${SAMPLES[@]}; do extract_length_fastq $i"_ccs.fastq" $i"_ccs" $CCS $LENGTHS; done
#for i in ${SAMPLES[@]}; do extract_length_fastq $i".clustered.hq.fastq" $i"_clustered" $INDIVIDUAL_CLUSTERED $LENGTHS; done
#for i in ${SAMPLES[@]}; do extract_length_fastq $i".collapsed.filtered.rep.fq" $i"_collapsed_filtered" $INDIVIDUAL_TOFU $LENGTHS; done

# merged samples 
#for i in {0..2}; do extract_length_fastq ${SAMPLES_MERGED[$i]}".clustered.hq.fastq" ${SAMPLES_MERGED[$i]}"_clustered" ${MERGED_CLUSTERED[$i]} $LENGTHS; done 
#for i in {0..2}; do extract_length_fastq ${SAMPLES_MERGED[$i]}".collapsed.filtered.rep.fq" ${SAMPLES_MERGED[$i]}"_collapsed_filtered" ${MERGED_TOFU[$i]} $LENGTHS; done 

# move output files to $LENGTHS 
#for i in { $INDIVIDUAL_CLUSTERED $INDIVIDUAL_TOFU $MERGED_CLUSTERED $MERGED_TOFU }; do cd $i; mv *seqlengths.txt* $LENGTHS; done

# convert raw subreads bam file 
for i in $RAW_DIR/*subreads.bam$*; do 
  sample=$(basename "$i" | cut -d "." -f 1,2 )
  echo $sample
  convert_bam2fq $sample "" $RAW_DIR $RAW_DIR
  #extract_length_fastq $sample".fastq" $sample"_subreads" $RAW_DIR $LENGTHS
done

cd $LENGTHS; mv *log* Log_output

