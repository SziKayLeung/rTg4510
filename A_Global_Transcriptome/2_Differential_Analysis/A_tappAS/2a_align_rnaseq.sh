#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=STAR_All-%A_%a.o
#SBATCH --array=0-63%8


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/2_Differential_Analysis/01_source_functions.sh


##-------------------------------------------------------------------------

SAMPLE=${RNASEQ_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}

# run_star $sample $Tg4510/J20/input_directory $MAPPED_output_directory  
run_star ${SAMPLE} ${RNASEQ_FILTERED_DIR} $RNASEQ_MAPPED_DIR/${sample} 