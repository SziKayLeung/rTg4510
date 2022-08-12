#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=ADBDR_Rnaseq_Star-%A_%a.o
#SBATCH --error=ADBDR_Rnaseq_Star-%A_%a.e
#SBATCH --array=0-63%8 # 64 samples

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

# run as array (defined in config file)
SAMPLE=${RNASEQ_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}


##-------------------------------------------------------------------------
echo "#************************************* RNAseq [Function 9, 10]"
## 9) run_star <list_of_samples> <rnaseq_input_directory> <output_dir>
run_star ${SAMPLE} $RNASEQ_FILTERED_DIR $RNASeq_MAPPED_DIR