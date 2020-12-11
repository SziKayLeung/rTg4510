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
#SBATCH --array=0-11 # 12 samples
#SBATCH --output=Isoseq3_alls-%A_%a.o
#SBATCH --error=Isoseq3_all-%A_%a.e

# 11-14/04/2019: run Isoseq3.1.2 for all Tg4510 samples: L22, K18, O23, S18, K17 
# 27/09/2019: run Isoseq3.2.2 for Samples 1-16 defined in raw.txt 
# 21/10/2020: Rerun with Isoseq3.4, ccs5.0, lima2.0 
# 03/11/2020: extract stats from CCS and LIMA

#************************************* DEFINE GLOBAL VARIABLES
# File directories 
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/Isoseq3_WKD
#cd $Isoseq3_WKD
#mkdir CCS LIMA REFINE CLUSTER
CCS=$Isoseq3_WKD/CCS
LIMA=$Isoseq3_WKD/LIMA
REFINE=$Isoseq3_WKD/REFINE
CLUSTER=$Isoseq3_WKD/CLUSTER
All_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
TG_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/TG_Merged
WT_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/WT_Merged

# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
SAMPLES_NAMES=(O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23)
SAMPLE=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
cd $FUNCTIONS
cat Isoseq_MouseRaw.txt
# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
BAM_FILES=(`cat "Isoseq_MouseRaw.txt" | egrep -v "^\s*(#|$)"`)
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}

#************************************* TO RUN FUNCTIONS ON WORKING SCRIPT
source $FUNCTIONS/Isoseq3.2.2_Functions.sh

#************************************* ALL samples separately 
# Isoseq3.4.0
    # run_CCS_batch <input_ccs_bam> <prefix_output_name> <Output_directory>
    # run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
    # run_REFINE $Sample $Input_LIMA_directory $Output_directory
    # run_CLUSTER $Sample $Input_REFINE_directory $Output_directory 
run_CCS_batch ${BAM_FILE} ${SAMPLE} $CCS
run_LIMA ${SAMPLE} $CCS $LIMA "no_multiplex"
run_REFINE ${SAMPLE} $LIMA $REFINE 
run_CLUSTER ${SAMPLE} $REFINE $CLUSTER

python $FUNCTIONS/Run_Stats/CCS.py $CCS "" All
python $FUNCTIONS/Run_Stats/LIMA.py $LIMA "" All

#************************************* Individual samples

# convert_fa2fq <file_name> <input_dir>
# run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
# tofu <prefix_sample> <cluster_dir> <input_dir> <output_dir>

convert_fa2fq ${SAMPLES}.clustered.hq.fasta $CLUSTER
run_minimap2 ${SAMPLES} $CLUSTER mm10 $Individual_Isoseq3_WKD/MAPPING
tofu ${SAMPLES} $CLUSTER $Individual_Isoseq3_WKD/MAPPING $Individual_Isoseq3_WKD/TOFU

