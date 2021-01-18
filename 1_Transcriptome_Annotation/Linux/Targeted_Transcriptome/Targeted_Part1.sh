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
#SBATCH --array=0-3 # 4 samples
#SBATCH --output=Targeted_Part1-%A_%a.o
#SBATCH --error=Targeted_Part1-%A_%a.e

# 11/12/2020: IsoSeq Analysis pipeline (Part1) for all targeted mouse batches (pooled and individual): IsoSeq3, Mapping, TOFU

#************************************* DEFINE VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq
Post_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq
cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER
cd $Post_Isoseq3_WKD; mkdir MAPPING TOFU CHAIN SQANTI
Targeted_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/

# For Pooled Targeted
### Important order of BAM files is the same as the sample names
SAMPLES_NAMES=(TargetedSeq1 Targeted_Seq_2 Targeted_Seq_3a Targeted_Seq_3b)
cat $Targeted_dir/Isoseq_Targeted_MouseRaw.txt
BAM_FILES=(`cat $Targeted_dir/Isoseq_Targeted_MouseRaw.txt | egrep -v "^\s*(#|$)"`)


#************************************* Pooled Targeted: Run IsoSeq Analysis pipeline
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
# Isoseq3.4.0
    # run_CCS_batch <input_ccs_bam> <prefix_output_name> <Output_directory>
    # run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
    # run_REFINE $Sample $Input_LIMA_directory $Output_directory
    # run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CCS_batch ${BAM_FILE} ${SAMPLE} $Isoseq3_WKD/CCS
run_LIMA ${SAMPLE} $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA multiplex
