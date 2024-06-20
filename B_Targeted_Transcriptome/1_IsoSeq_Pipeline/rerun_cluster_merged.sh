#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Targeted_Mouse_Part3.o
#SBATCH --error=Targeted_Mouse_Part3.e

# 03/06/2021: Same script as Targeted_Mouse_Part2.sh but only work with the same samples as Whole Transcriptome for fair comparisons
# 14/01/2022: Re-run SQANIT3 on subset of whole samples with no RNA-Seq Support
##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_isoseq_processing

export samplename=MatchedMouse
mkdir -p $WKD_ROOT/7b_matched_only

cd $WKD_ROOT/7b_matched_only/8_merged_cluster
source activate isoseq3
isoseq3 cluster $samplename.flnc.fofn $samplename.clustered.bam --verbose --use-qvs