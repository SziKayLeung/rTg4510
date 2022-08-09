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
#SBATCH --reservation=research_project-mrc148213_5
#SBATCH --output=Targeted_Mouse_Part1.o
#SBATCH --error=Targeted_Mouse_Part1.e

# 07/04/2021: Inital script outline

#************************************* DEFINE GLOBAL VARIABLES
# File directories
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse; mkdir IsoSeq Post_IsoSeq RNASeq
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/RNASeq

cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER
cd $Isoseq3_WKD/REFINE; mkdir BATCHES
cd $Isoseq3_WKD/CLUSTER; mkdir BATCHES
cd $PostIsoseq3_WKD; mkdir mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER
cd $RNASeq_WKD; mkdir MAPPED

# For Pooled Targeted
### Important order of BAM files is the same as the sample names
BATCH_NAMES=(Targeted_Seq_1 Targeted_Seq_2 Targeted_Seq_3a Targeted_Seq_3b)
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome
cat $RAWDIR/Isoseq_Targeted_MouseRaw.txt
BAM_FILES=(`cat $RAWDIR/Isoseq_Targeted_MouseRaw.txt | egrep -v "^\s*(#|$)"`)

# Other input files and directory
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
DEMUX_SCRIPT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome/

# sourcing functions script
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome
source $FUNCTIONS/Targeted_Mouse_Functions.sh

module load Miniconda2/4.3.21
################################################################################################
echo "#*************************************  Isoseq3 [Function 1, 2]"
BATCH=${BATCH_NAMES[${SLURM_ARRAY_TASK_ID}]}
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}
# Isoseq3.4.0
    # run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory>
    # run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
    # run_REFINE $Sample $Input_LIMA_directory $Output_directory
    # run_targeted_REFINE $Input_Pooled_Sample $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
    # run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CCS ${BAM_FILE} ${BATCH} $Isoseq3_WKD/CCS
run_LIMA ${BATCH} $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA "no_multiplex"
run_REFINE ${BATCH} $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE/BATCHES
run_CLUSTER ${BATCH} $Isoseq3_WKD/REFINE/BATCHES $Isoseq3_WKD/CLUSTER/BATCHES
