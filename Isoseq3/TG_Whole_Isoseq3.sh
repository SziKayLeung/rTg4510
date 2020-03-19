#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 12/02/2020: run merged Isoseq3.2.2 on TG6 samples (Tg4510; 3 samples at 2mos and 3 samples at 8mos)
    # Samples O18, K18, S18, L22, Q20, K24

#************************************* DEFINE GLOBAL VARIABLES
sl693=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693
WKD=$sl693/WholeTranscriptome/Tg4510/WT_vs_TG/TG_all
Isoseq3_WKD=$WKD/Isoseq3_WKD
RAW_DIR=$sl693/WholeTranscriptome/rawdata


#************************************* TO RUN FUNCTIONS ON WORKING SCRIPT
# sourcing functions script and input directories
module load Miniconda2
source $sl693/Scripts/general/IsoSeq/Isoseq3.2.2_Functions.sh
source activate isoseq3

## If working on merging specific samples 
# DEFINE SAMPLES_NAMES list in working script
# DEFINE RAW_DATA_SUBREADS in working script (path directories of subreads.bam files)
# run_isoseq3.2.1_merge <Sample_name_output> <WKD>

SAMPLES_NAMES=(O18 K18 S18 L22 Q20 K24)

RAW_DATA_SUBREADS=(
#2.O18_Tg4510_TG_2mos
$RAW_DIR/m54082_180605_141944.subreads.bam
#8.K18_Tg4510_TG_2months
$RAW_DIR/m54082_190307_045507.subreads.bam
#11.S18_Tg4510_TG_2mos
$RAW_DIR/m54082_190404_101400.subreads.bam
#7.L22_Tg4510_TG_8months
$RAW_DIR/m54082_190306_083150.subreads.bam
#15.Q20_Tg4510_TG_8mos 
$RAW_DIR/m54082_190527_173356.subreads.bam
#16.K24_Tg4510_TG_8mos
$RAW_DIR/m54082_190529_082942.subreads.bam
)

run_isoseq3.2.1_merge TG_Whole $Isoseq3_WKD
