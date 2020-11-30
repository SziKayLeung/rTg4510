#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 12/02/2020: run merged Isoseq3.2.2 on WT6 samples (Tg4510; 3 samples at 2mos and 3 samples at 8mos)
    # Samples Q21, K17, M21, O23, S23, K23 

#************************************* DEFINE GLOBAL VARIABLES
sl693=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693
WKD=$sl693/WholeTranscriptome/Tg4510/WT_vs_TG/WT_all
Isoseq3_WKD=$WKD/Isoseq3_WKD
RAW_DIR=$sl693/WholeTranscriptome/rawdata


#************************************* TO RUN FUNCTIONS ON WORKING SCRIPT
# sourcing functions script and input directories
source $sl693/Scripts/general/IsoSeq/Isoseq3.2.2_Functions.sh

## If working on merging specific samples 
# DEFINE SAMPLES_NAMES list in working script
# DEFINE RAW_DATA_SUBREADS in working script (path directories of subreadset.xml files)
# run_isoseq3.2.1_merge <Sample_name_output> <WKD>

SAMPLES_NAMES=(Q21 K17 M21 O23 S23 K23)

RAW_DATA_SUBREADS=(
#1. Q21_Tg4510_WT_2mos
$RAW_DIR/m54082_180607_173058.subreads.bam
#4. K17_Tg4510_WT_2mos
$RAW_DIR/m54082_190405_063832.subreads.bam
#11. M21_Tg4510_WT_2mos
$RAW_DIR/m54082_190430_163756.subreads.bam
#2. O23_Tg4510_WT_8mos
$RAW_DIR/m54082_190401_165425.subreads.bam
#3. S23_Tg4510_WT_8mos
$RAW_DIR/m54082_190403_135102.subreads.bam
#12. K23_Tg4510_WT_8mos
$RAW_DIR/m54082_190524_145911.subreads.bam
)

run_isoseq3.2.1_merge WT_Whole $Isoseq3_WKD
