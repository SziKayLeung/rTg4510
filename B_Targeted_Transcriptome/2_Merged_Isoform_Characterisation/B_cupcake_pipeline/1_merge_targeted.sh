#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../../../bash_output/1_merge_targeted.o
#SBATCH --error=../../../bash_output/1_merge_targeted.e

##-------------------------------------------------------------------------

SC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/B_cupcake_pipeline

# first job - no dependencies 
# copy samples to folder and rename
jid1=$(sbatch ${SC}/1a_prepare_samples.sh)

# second job - run after jid finishes successfully
# batch and align samples separately
jid2=$(sbatch  --dependency=afterok:$jid1 ${SC}/1b_batch_align_filter.sh)

# third job - run after jid2 finishes successfully
# merge all samples, collapse and run sqanti
jid3=$(sbatch --dependency=afterany:$jid2 --job-name=collapse ${SC}/1c_merged_collapse_sqanti3.sh)