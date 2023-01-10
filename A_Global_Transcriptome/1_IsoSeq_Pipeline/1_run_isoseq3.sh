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
#SBATCH --array=0-11 # 12 samples
#SBATCH --output=All_Batch_Processing-%A_%a.o
#SBATCH --error=All_Batch_Processing-%A_%a.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

# run as array (defined in config file)
SAMPLE=${ALL_SAMPLE_NAMES[${SLURM_ARRAY_TASK_ID}]}
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}


##-------------------------------------------------------------------------

# Isoseq3.4.0
# run_CCS_batch <input_ccs_bam> <prefix_output_name> <Output_directory>
# run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
# run_REFINE $Sample $Input_LIMA_directory $Output_directory
# run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CCS ${BAM_FILE} ${SAMPLE} ${WKD_ROOT}/1_isoseq3/1_ccs
run_LIMA ${SAMPLE} ${WKD_ROOT}/1_isoseq3/1_ccs ${WKD_ROOT}/1_isoseq3/2_lima "no_multiplex"
run_REFINE ${SAMPLE} ${WKD_ROOT}/1_isoseq3/2_lima ${WKD_ROOT}/1_isoseq3/3_refine
run_CLUSTER ${SAMPLE} ${WKD_ROOT}/1_isoseq3/3_refine ${WKD_ROOT}/1_isoseq3/4_cluster


##-------------------------------------------------------------------------
run_star ${SAMPLE} ${RNASEQ_FILTERED_DIR} ${RNASEQ_MAPPED_DIR}