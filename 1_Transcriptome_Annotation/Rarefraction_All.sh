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
#SBATCH --array=0-2
#SBATCH --output=rarefaction.o
#SBATCH --error=rarefaction.e

module load Miniconda2 
source activate sqanti2_py3

#************************************* DEFINE GLOBAL VARIABLES
### IMPORTANT TO ENSURE ORDER OF SAMPLES AND DIRECTORY: ALL_Merged, TG_Merged, WT_Merged ###
SAMPLES=(All_Merged TG_Merged WT_Merged)
DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/
TOFU_INPUT_DIR=($DIR/All_Merged/TOFU/mm10 $DIR/TG_Merged/TOFU/mm10 $DIR/WT_Merged/TOFU/mm10)
RAREFACTION_OUTPUT_DIR=($DIR/All_Merged/RAREFACTION $DIR/TG_Merged/RAREFACTION $DIR/WT_Merged/RAREFACTION)
SQANTI_INPUT_DIR=($DIR/All_Merged/SQANTI $DIR/TG_Merged/SQANTI $DIR/WT_Merged/SQANTI)

source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh
collapsed.filtered_classification.filtered_lite_classification.txt

# submit to 3 batches with a counting track looping in through each sample, input and output directory
# make_file_for_rarefaction <sample_name_prefix> <input_tofu_directory> <working_directory> <input_sqanti_directory>
# note: input_file = <sample_name_prefix>.collapsed.filtered; run after tofu and SQANTI

count=(1 2 3)
count_track=${count[${SLURM_ARRAY_TASK_ID}]}
make_file_for_rarefaction ${SAMPLES[$count_track]} ${TOFU_INPUT_DIR[$count_track]} ${RAREFACTION_OUTPUT_DIR[$count_track]} ${SQANTI_INPUT_DIR[$count_track]]}; done 
