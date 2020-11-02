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
#SBATCH --output=Post_Isoseq3_all.output

# 22/10/2020: Merging All, TG and WT samples after Isoseq3_all.sh

#************************************* DEFINE GLOBAL VARIABLES
# File directories 
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/Isoseq3_WKD
REFINE=$Isoseq3_WKD/REFINE
All_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
TG_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/TG_Merged
WT_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/WT_Merged

#************************************* TO RUN FUNCTIONS ON WORKING SCRIPT
#source $FUNCTIONS/Isoseq3.2.2_Functions.sh
merging_at_refine(){
  module load Miniconda2/4.3.21
  source activate isoseq3
  isoseq3 --version #isoseq3 3.2.2 (commit v3.2.2)

  ###********************* Merging at REFINE 
  # Define variable "Merge_Samples" as a list of all samples, in order to find the specified flnc.bam (for dataset create ConsensusReadSet) 
  # Define variable "all_flnc_bams" for merged path-directory of all flnc samples (for dataset create ConsensusReadSet)   
  Merge_Samples=$(echo "${@:4}")

  echo "Merging flnc of samples $Merge_Samples"
  all_flnc_bams=$(
      for i in ${Merge_Samples[@]}; do
          flnc_bam_name=$(find $1 -name "*.flnc.bam" -exec basename \{} \; | grep ^$i )
          flnc_bam=$(find $1 -name "*.flnc.bam" | grep "$flnc_bam_name" )
          echo $flnc_bam
      done
  )
  
  cd $2
  printf '%s\n' "${all_flnc_bams[@]}" > $3.flnc.fofn
  cat $.flnc.fofn
  
  ###*********************
  
  isoseq3 cluster $3.flnc.fofn $3.clustered.bam --verbose --use-qvs
  gunzip *.gz*
  # convert clustered.hq.fasta to clustered.hq.fastq for 
  source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh
  convert_fa2fq $3.clustered.hq.fasta .

  source deactivate
}

#************************************* All samples, WT, TG merged at refine  
# merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
merging_at_refine $REFINE $All_Isoseq3_WKD/CLUSTER All_Merged O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23
merging_at_refine $REFINE $TG_Isoseq3_WKD/CLUSTER TG_Merged O18 K18 S18 L22 Q20 K24
merging_at_refine $REFINE $WT_Isoseq3_WKD/CLUSTER WT_Merged Q21 K17 M21 O23 S23 K23