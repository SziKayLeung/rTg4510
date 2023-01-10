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
#SBATCH --array=0-2 # 3 samples
#SBATCH --output=Targeted_ADBDR_Part1-%A_%a.o
#SBATCH --error=Targeted_ADBDR_Part1-%A_%a.e

# 07/04/2021: Inital script outline

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

# run as array (defined in config file)
BATCH=${BATCH_NAMES[${SLURM_ARRAY_TASK_ID}]}
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}

##-------------------------------------------------------------------------
echo "#*************************************  Isoseq3 [Function 1, 2]"
# Isoseq3.4.0
    # run_CCS <sample>
    # run_LIMA <sample> <"no_multiplex"/"multiplex">
    # run_REFINE <sample> 
    # run_CLUSTER <sample>
run_CCS ${BAM_FILE} ${BATCH} 
run_LIMA ${BATCH} "multiplex"
run_LIMA ${BATCH} "no_multiplex"
run_REFINE ${BATCH} 
run_CLUSTER ${BATCH} 