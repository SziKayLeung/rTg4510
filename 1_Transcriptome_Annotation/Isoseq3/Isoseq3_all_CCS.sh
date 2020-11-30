#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=20:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 22/04/2020: re_run Isoseq3.2.2 for Samples 1-16 defined in raw.txt to generate CCS output.txt

#************************************* DEFINE GLOBAL VARIABLES
# File directories 
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD
#cd $Isoseq3_WKD
#mkdir CCS LIMA REFINE CLUSTER
CCS=$Isoseq3_WKD/CCS
LIMA=$Isoseq3_WKD/LIMA
REFINE=$Isoseq3_WKD/REFINE
CLUSTER=$Isoseq3_WKD/CLUSTER

# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
SAMPLES_NAMES=(Q21)
cd $FUNCTIONS
cat Isoseq_MouseRaw.txt
# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
BAM_FILES=(`cat "Isoseq_MouseRaw.txt" | egrep -v "^\s*(#|$)"`)

#************************************* TO RUN FUNCTIONS ON WORKING SCRIPT
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
for i in ${SAMPLES_NAMES[@]}; do
    echo $i
    run_CCS $CCS
    count=$((count+1))
done 
