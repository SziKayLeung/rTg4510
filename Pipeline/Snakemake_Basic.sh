#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=4:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=SnakeMakeCheck_rerun.o
#SBATCH --output=SnakeMakeCheck_rerun.e

# 11/12/2020: Run 2 samples to test output of IsoSeq and Post IsoSeq Pipeline

module load Miniconda2
source activate sqanti2_py3
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/Pipeline/
snakemake -j 16 --use-conda -s Snakefile_Basic --configfile config_basic.yaml
