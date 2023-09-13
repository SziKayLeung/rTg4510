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
#SBATCH --output=../../../bash_output/2_post_align_rnaseq2_isoseq.o
#SBATCH --error=../../../bash_output/2_post_align_rnaseq2_isoseq.e

##-------------------------------------------------------------------------

SC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome/2_Differential_Analysis/B_DESeq2

# first job - no dependencies 
# align rnaseq individual samples to isoseq collapsed annotation
echo "Aligning RNA-Seq files to Iso-Seq defined transcriptome using kallisto"
jid1=$(sbatch ${SC}/2a_align_rnaseq2isoseq.sh)

# second job - run after jid finishes successfully
# merge alignment reads
jid2=$(sbatch  --dependency=afterany:$jid1 ${SC}/2b_merge_rnaseq_alignment.sh)