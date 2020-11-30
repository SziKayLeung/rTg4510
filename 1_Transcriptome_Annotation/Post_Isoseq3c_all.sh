#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-11 # 12 samples

# 27/10/2020: Minimap2, TOFU Merged Inidivudal IsoSeq samples

#************************************* DEFINE GLOBAL VARIABLES
Individual_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/Isoseq3_WKD/CLUSTER

### setting array names for samples
SAMPLES_BATCHES=(O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23) 
SAMPLES=${SAMPLES_BATCHES[${SLURM_ARRAY_TASK_ID}]}

source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh

#************************************* Individual samples

# convert_fa2fq <file_name> <input_dir>
# run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
# tofu <prefix_sample> <cluster_dir> <input_dir> <output_dir>

convert_fa2fq ${SAMPLES}.clustered.hq.fasta $CLUSTER
run_minimap2 ${SAMPLES} $CLUSTER mm10 $Individual_Isoseq3_WKD/MAPPING
tofu ${SAMPLES} $CLUSTER $Individual_Isoseq3_WKD/MAPPING $Individual_Isoseq3_WKD/TOFU
