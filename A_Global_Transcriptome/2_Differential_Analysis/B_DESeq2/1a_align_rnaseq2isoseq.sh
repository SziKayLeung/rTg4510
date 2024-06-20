#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-58%8 #59 samples (64-5)
#SBATCH --output=../../../bash_output/2b_align-rnaseq2isoseq-%A_%a.o
#SBATCH --error=../../../bash_output/2b_align-rnaseq2isoseq-%A_%a.e

## ---------------------------
## Purpose: Align RNA-Seq to Iso-Seq annotation for differential expression analysis
## 
## 24/02/2023: Align RNA-Seq samples to Iso-Seq defined transcriptome
## ---------------------------

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/2_Differential_Analysis/01_source_functions.sh
source $SC_ROOT/2_Differential_Analysis/rTg4510_differential.config

RNASEQ_SAMPLES_NAMES=(${RNASEQ_SAMPLES_NAMES[@]})
SAMPLE=${RNASEQ_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}


##-------------------------------------------------------------------------

run_kallisto_1sample ${RNASEQ_FILTERED_DIR} ${SAMPLE} ${NAME}_Kallisto.idx $WKD_ROOT/2_post_isoseq3/10_rnaseq2isoseq