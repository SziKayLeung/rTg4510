#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=00:10:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

# aim: generate a file of the missing ONT reads that are documented in the gtf but not in the abundance file
# motivation:
  # noticed some transcripts in TALON gtf but not in the TALON abundance file. 
  # output these transcripts

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/rTg4510_merged.config
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/01_source_function.sh


##-------------------------------------------------------------------------

source activate sqanti2_py3

# difference between talon abundance and gtf 
Rscript ${LOGEN_ROOT}/assist_ont_processing/identify_talon_reads_no_counts.R \
  -g ${ONT_UNFILTERED_DIR}/ONTTargeted_unfiltered_talon.gtf \
  -c ${ONT_UNFILTERED_DIR}/ONTTargeted_unfiltered_talon_abundance.tsv \
  -o ${ONT_UNFILTERED_DIR}/ONTTargeted_unfiltered_talon_missing_abundance.txt