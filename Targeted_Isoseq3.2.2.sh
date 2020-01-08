#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrcq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 08/01/2019: Created script to run Isoseq3.2.2 on Mouse Targeted Run 1 (TargetedSeq1)
    # Samples K20, K18, K23, K21, K19, K17 
    # Aaron transferred *subreads* files 

sl693=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693
FUNCTIONS=$sl693/Scripts/general/IsoSeq

#************************************* DOWNLOAD AND PREPARE RAWDATA
# Prepare subread file 
cd $FUNCTIONS
cat >> Isoseq_Targeted_MouseRaw.txt <<EOL
#1.Targeted_Sequencing_1: K20, K18, K23, K21, K19, K17 
/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/rawdata/m54082_191116_131337.subreads.bam
EOL

cat >> Isoseq_Targeted_rawdir.txt <<EOL
#1.Targeted_Sequencing_1: K20, K18, K23, K21, K19, K17 
/bioseqfs/PacBio/Sequel/3087_10040/r54082_20191115_164308/2_B01/m54082_191116_131337.subreadset.xml
/bioseqfs/PacBio/Sequel/3087_10040/r54082_20191115_164308/2_B01/m54082_191116_131337.subreads.bam
/bioseqfs/PacBio/Sequel/3087_10040/r54082_20191115_164308/2_B01/m54082_191116_131337.subreads.bam.pbi
EOL


#************************************* DEFINE VARIABLES
TS1_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/Targeted_Seq_1
SAMPLES_NAMES=(TargetedSeq1)
BAM_FILES=(`cat $FUNCTIONS/Isoseq_Targeted_MouseRaw.txt | egrep -v "^\s*(#|$)"`)


#************************************* Run script
git checkout remotes/origin/Isoseq3.2.2 
source $FUNCTIONS/Isoseq3.2.2_Functions.sh

run_CCS $TS1_dir/CCS
for i in ${SAMPLES_NAMES[@]}; do  
    run_LIMA $i $TS1_dir/CCS $TS1_dir/LIMA multiplex
    run_REFINE $i $TS1_dir/LIMA $TS1_dir/REFINE 
    run_CLUSTER $i $TS1_dir/REFINE $TS1_dir/CLUSTER
done
