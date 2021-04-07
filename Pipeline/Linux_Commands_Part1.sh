#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Linux_Commands-%A_%a.o
#SBATCH --error=Linux_Command-%A_%a.e
#SBATCH --array=0-1

# 11/12/2020: Run 2 samples to test output of IsoSeq and Post IsoSeq Pipeline

#************************************* DEFINE GLOBAL VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Testing/Linux
cd $Isoseq3_WKD
mkdir CCS LIMA REFINE CLUSTER MAPPING TOFU
CCS=$Isoseq3_WKD/CCS
LIMA=$Isoseq3_WKD/LIMA
REFINE=$Isoseq3_WKD/REFINE
CLUSTER=$Isoseq3_WKD/CLUSTER
MAPPING=$Isoseq3_WKD/MAPPING
TOFU=$Isoseq3_WKD/TOFU

# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
SAMPLES_NAMES=(O18 S18)
SAMPLE=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
cat Isoseq_MouseRaw_Testing.txt
# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
BAM_FILES=(`cat "Isoseq_MouseRaw_Testing.txt" | egrep -v "^\s*(#|$)"`)
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}

#************************************* IsoSeq - ALL samples separately
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
# Isoseq3.4.0
    # run_CCS_batch <input_ccs_bam> <prefix_output_name> <Output_directory>
    # run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
    # run_REFINE $Sample $Input_LIMA_directory $Output_directory
    # run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CCS_batch ${BAM_FILE} ${SAMPLE} $CCS
run_LIMA ${SAMPLE} $CCS $LIMA "no_multiplex"
run_REFINE ${SAMPLE} $LIMA $REFINE
run_CLUSTER ${SAMPLE} $REFINE $CLUSTER

#************************************* Post-IsoSeq - ALL samples separately
source $FUNCTIONS/Post_Isoseq3_Function.sh
# convert_fa2fq <file_name> <input_dir>
# run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
# tofu <prefix_sample> <cluster_dir> <input_dir> <output_dir>
convert_fa2fq ${SAMPLE}.clustered.hq.fasta $CLUSTER
run_minimap2 ${SAMPLE} $CLUSTER mm10 $MAPPING
tofu ${SAMPLE} $CLUSTER $MAPPING $TOFU
