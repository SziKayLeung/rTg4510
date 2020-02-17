#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 06/10/2019: Created Script to run RNASeq, FeatureCounts and Post_Isoseq3_Functions on Samples 1-16
    # note all STAR and featurecounts mapping done with gencode.vM22.annotation.gtf
# 17/02/2020: Reran SQANTI2 (v7.2)

#************************************* DEFINE GLOBAL VARIABLES
# File directories 
POLISHED=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER # same folder as CLUSTER, used as reference for functions script to contain hq.fasta
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
FEATURECOUNTS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/FeatureCounts/Whole_Transcriptome
MAPPING=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/MAPPING
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/TOFU
SQANTI2_output_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/SQANTI2
STAR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/MAPPED/Individual

#************************************* TO RUN FUNCTIONS ON WORKING SCRIPT

module load Miniconda2/4.3.21
source activate sqanti2_py3

SAMPLES_NAMES=(Q21 O18 C21 E18 C20 B21 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)

# RNASeq
#source $FUNCTIONS/RNASeq/STAR_Functions.sh
#for i in ${SAMPLES_NAMES[@]}; do
    # run_star $sample $Tg4510/J20/input_directory $MAPPED_output_directory  $REFERENCE 
    #run_star $i $RNASeq_Filtered $STAR $REFERENCE
#done 

# FeatureCounts
# run_featurecounts at TRANSCRIPT level of all specified samples
# run_featurecounts_transcript_specified <input_dir> <input_reference_dir> <output_prefix_name> <output_dir>
# <input_dir> containing mapped, sorted bam files from STAR (RNASeq)
#source $FUNCTIONS/RNASeqvsIsoseq/Run_FeatureCounts.sh
#run_featurecounts_transcript_specified $STAR $REFERENCE All_Whole_Transcriptome $FEATURECOUNTS

# Post_IsoSeq3
#cd $POLISHED; gunzip *.gz
source $FUNCTIONS/Post_IsoSeq/Post_Isoseq3_Functions.sh
for i in ${SAMPLES_NAMES[@]}; do 
    #convert_fa2fq $i $POLISHED
    #run_minimap2 $i
    #tofu $i
    run_sqanti2_QC $i
    run_sqanti2_Filter $i 
done

source deactivate
