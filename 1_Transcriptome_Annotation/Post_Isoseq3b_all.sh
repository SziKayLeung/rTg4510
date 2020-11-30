#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-2 # 3 samples
#SBATCH --output=Post_IsoSeq3b-%A_%a.o
#SBATCH --error=Post_IsoSeq3b-%A_%a.e

# 23/10/2020: Minimap2, TOFU Merged, WT vs TG samples 

#************************************* DEFINE GLOBAL VARIABLES
All_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
TG_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/TG_Merged
WT_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/WT_Merged

source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh

#cd $All_Isoseq3_WKD; mkdir MAPPING TOFU; cd MAPPING; mkdir ERCC mm10; cd .. ; cd TOFU; mkdir ERCC mm10 
#cd $WT_Isoseq3_WKD; mkdir MAPPING TOFU; cd MAPPING; mkdir ERCC mm10; cd .. ; cd TOFU; mkdir ERCC mm10 
#cd $TG_Isoseq3_WKD; mkdir MAPPING TOFU; cd MAPPING; mkdir ERCC mm10; cd .. ; cd TOFU; mkdir ERCC mm10 

### setting array names for samples
SAMPLES_BATCHES=(All_Merged TG_Merged WT_Merged) 
SAMPLES_BATCHES_CLUSTER_DIR=($All_Isoseq3_WKD/CLUSTER $TG_Isoseq3_WKD/CLUSTER $WT_Isoseq3_WKD/CLUSTER)
SAMPLES_BATCHES_MM10_MAPPING_DIR=($All_Isoseq3_WKD/MAPPING/mm10 $TG_Isoseq3_WKD/MAPPING/mm10 $WT_Isoseq3_WKD/MAPPING/mm10)
SAMPLES_BATCHES_ERCC_MAPPING_DIR=($All_Isoseq3_WKD/MAPPING/ERCC $TG_Isoseq3_WKD/MAPPING/ERCC $WT_Isoseq3_WKD/MAPPING/ERCC)
SAMPLES_BATCHES_MM10_TOFU_DIR=($All_Isoseq3_WKD/TOFU/mm10 $TG_Isoseq3_WKD/TOFU/mm10 $WT_Isoseq3_WKD/TOFU/mm10)
SAMPLES_BATCHES_ERCC_TOFU_DIR=($All_Isoseq3_WKD/TOFU/ERCC $TG_Isoseq3_WKD/TOFU/ERCC $WT_Isoseq3_WKD/TOFU/ERCC)

SAMPLES=${SAMPLES_BATCHES[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_CLUSTER_DIR=${SAMPLES_BATCHES_CLUSTER_DIR[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_MM10_MAPPING_DIR=${SAMPLES_BATCHES_MM10_MAPPING_DIR[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_ERCC_MAPPING_DIR=${SAMPLES_BATCHES_ERCC_MAPPING_DIR[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_MM10_TOFU_DIR=${SAMPLES_BATCHES_MM10_TOFU_DIR[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_ERCC_TOFU_DIR=${SAMPLES_BATCHES_ERCC_TOFU_DIR[${SLURM_ARRAY_TASK_ID}]}

#************************************* All samples, WT, TG Minimap 

# convert_fa2fq <file_name> <input_dir>
convert_fa2fq All_Merged.clustered.hq.fasta $All_Isoseq3_WKD/CLUSTER
convert_fa2fq WT_Merged.clustered.hq.fasta $WT_Isoseq3_WKD/CLUSTER
convert_fa2fq TG_Merged.clustered.hq.fasta $TG_Isoseq3_WKD/CLUSTER

# run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
#run_minimap2 ${SAMPLES} ${SAMPLES_CLUSTER_DIR} mm10 ${SAMPLES_MM10_MAPPING_DIR} 
#run_minimap2 ${SAMPLES} ${SAMPLES_CLUSTER_DIR} ERCC ${SAMPLES_ERCC_MAPPING_DIR} 

#************************************* All samples, WT, TG TOFU

 # tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir>
tofu ${SAMPLES} ${SAMPLES_CLUSTER_DIR} ${SAMPLES_MM10_MAPPING_DIR} ${SAMPLES_MM10_TOFU_DIR}
tofu ${SAMPLES} ${SAMPLES_CLUSTER_DIR} ${SAMPLES_ERCC_MAPPING_DIR} ${SAMPLES_ERCC_TOFU_DIR}
