#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=10:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 17/10/2019: Created script to chain all Tg4510 samples according to cupcake 
## note modified chain_samples.py to take the same path but different sample name files for group_filename, gff_filename, count_filename, and fastq_filename 

#************************************* DEFINE GLOBAL VARIABLES
module load Miniconda2
source activate sqanti2

FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Modified_Scripts
CHAIN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/CHAIN
TOFU_PATH=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/TOFU

# Prepare configuration file 
cd $CHAIN
cat >> config_file <<EOL
SAMPLE=Q21;$TOFU_PATH
SAMPLE=O18;$TOFU_PATH
SAMPLE=C21;$TOFU_PATH
SAMPLE=E18;$TOFU_PATH
SAMPLE=C20;$TOFU_PATH
SAMPLE=B21;$TOFU_PATH
SAMPLE=L22;$TOFU_PATH
SAMPLE=K18;$TOFU_PATH
SAMPLE=O23;$TOFU_PATH
SAMPLE=S23;$TOFU_PATH
SAMPLE=S18;$TOFU_PATH
SAMPLE=K17;$TOFU_PATH
SAMPLE=M21;$TOFU_PATH
SAMPLE=K23;$TOFU_PATH
SAMPLE=Q20;$TOFU_PATH
SAMPLE=K24;$TOFU_PATH

GROUP_FILENAME=.collapsed.group.txt
GFF_FILENAME=.collapsed.filtered.gff
COUNT_FILENAME=.collapsed.filtered.abundance.txt
FASTQ_FILENAME=.collapsed.filtered.rep.fq
EOL

#************************************* DEFINE GLOBAL VARIABLES
python $FUNCTIONS/chain_samples.py $CHAIN/config_file count_fl 