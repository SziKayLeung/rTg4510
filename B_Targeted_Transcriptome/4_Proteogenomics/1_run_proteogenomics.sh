#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrchq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=1_run_proteogenomics.o
#SBATCH --error=1_run_proteogenomics.e

# 04/01/2024: run proteogenomics pipeline on rTg4510 combined targeted dataset

#-----------------------------------------------------------------------#
## print start date and time
echo Job started on:
date -u

module load Miniconda2
source activate nanopore
source /lustre/projects/Research_Project-MRC148213/sl693/scripts/LOGen/proteomics/proteogenomics.sh
source /lustre/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/4_Proteogenomics/rTg4510_proteomics.config

echo "#************************************* Collate and prepare long-read data"
collate_longread_processed
prepare_reference_tables
summarise_longread_data

echo "#************************************* Call open reading frames and classify proteins"
call_orf
determine_best_orf
refine_calledorf
classify_protein

echo "#***************All done!****************#"