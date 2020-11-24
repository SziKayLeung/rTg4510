#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrcq # submit to the serial queue
#PBS -l walltime=50:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 27/05/2020: Created script to run Isoseq3.1.2 on Mouse Targeted Run 1 (TargetedSeq1) for Isoseq3 paper 
    # Samples K20, K18, K23, K21, K19, K17 

sl693=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693
FUNCTIONS=$sl693/Scripts/general/IsoSeq

#************************************* DOWNLOAD AND PREPARE RAWDATA
#cd $FUNCTIONS
#cat >> Isoseq_Targeted_MouseRawSub.txt <<EOL
#1.Targeted_Sequencing_1: K20, K18, K23, K21, K19, K17 
#/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/rawdata/m54082_191116_131337.subreadset.xml
#EOL

#************************************* DEFINE VARIABLES
TS1_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/Targeted_Seq_1/Isoseq3_Paper
# cd $TS1_dir; mkdir CCS LIMA REFINE CLUSTER POLISH

SAMPLES_NAMES=(TargetedSeq1)
Pooled_Samples=(K19 K23 K21 K18 K20 K17)
BAM_FILES=(`cat $FUNCTIONS/Isoseq_Targeted_MouseRaw.txt | egrep -v "^\s*(#|$)"`)
SUB_FILES=(`cat $FUNCTIONS/Isoseq_Targeted_MouseRawSub.txt | egrep -v "^\s*(#|$)"`)

#************************************* Run script
source $FUNCTIONS/Isoseq3.1.2_Functions.sh

# Isoseq3 workflow for pooled samples in pooled dataset
#run_CCS $TS1_dir/CCS
count=0
for i in ${SAMPLES_NAMES[@]}; do  
  #run_LIMA $i $TS1_dir/CCS $TS1_dir/LIMA "no_multiplex"
  #run_REFINE $i $TS1_dir/LIMA $TS1_dir/REFINE 
  #run_CLUSTER $i $TS1_dir/REFINE $TS1_dir/CLUSTER
  run_POLISH $i $TS1_dir/CLUSTER $FUNCTIONS $TS1_dir/POLISH
  count=$((count+1))
done
