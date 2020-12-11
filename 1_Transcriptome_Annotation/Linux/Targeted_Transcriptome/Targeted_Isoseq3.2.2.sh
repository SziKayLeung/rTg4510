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
#SBATCH --output=Post_IsoSeq3b-%A_%a.o
#SBATCH --error=Post_IsoSeq3b-%A_%a.e

# 08/01/2019: Created script to run Isoseq3.2.2 on Mouse Targeted Run 1 (TargetedSeq1)
    # Samples K20, K18, K23, K21, K19, K17 
    # Aaron transferred *subreads* files 
# 16/01/2019: Isoseq3 workflow for individual samples in pooled dataset
# 03/11/2020: Isoseq3 workflow for all targete mouse batches (pooled)

Rawdata=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/targeted_isoseq/Raw_Data

#************************************* DEFINE VARIABLES
TS1_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/Targeted_Seq_1
TS2_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/Targeted_Seq_2
TS3a_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/Targeted_Seq_3a
TS3b_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/Targeted_Seq_3b
TS_output_dir=($TS1_dir $TS2_dir $TS3a_dir $TS3b_dir)
SAMPLES_NAMES=(TargetedSeq1 Targeted_Seq_2 Targeted_Seq_3a Targeted_Seq_3b)

Pooled_Samples_Targeted1=(K19 K23 K21 K18 K20 K17)
Pooled_Samples_Targeted2=(S19 K24 L22 M21 O18 O23 O22 P19 T20)
Pooled_Samples_Targeted3a=(Q20 Q21 S18 S23 Q18 Q17 L18 Q23 T18) # failed 3rd run
Pooled_Samples_Targeted3b=(Q20 Q21 S18 S23 Q18 Q17 L18 Q23 T18) # successful 3rd run
BAM_FILES=(`cat $Rawdata/Isoseq_Targeted_MouseRaw.txt | egrep -v "^\s*(#|$)"`)
Barcoded_Targeted1_config_file=$Rawdata/Isoseq_Mouse_Targeted1_barcode.config
Barcoded_Targeted2_config_file=$Rawdata/Isoseq_Mouse_Targeted2_barcode.config
Barcoded_Targeted3_config_file=$Rawdata/Isoseq_Mouse_Targeted3_barcode.config

SAMPLES=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}
TS_output=${TS_output_dir[${SLURM_ARRAY_TASK_ID}]}

#************************************* Run script
#git checkout remotes/origin/Isoseq3.2.2 
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq/Isoseq3.2.2_Functions.sh

# create directories 
#cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/; mkdir Targeted_Seq_1 Targeted_Seq_2 Targeted_Seq_3a Targeted_Seq_3b
#for i in ${TS_output_dir[@]}; do cd $i; mkdir CCS LIMA REFINE CLUSTER; done

# Isoseq3 workflow for pooled samples in pooled dataset
# run_CCS_batch <input_ccs_bam> <prefix_output_name> <Output_directory>
# run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
run_CCS_batch ${BAM_FILE} ${SAMPLES} ${TS_output}/CCS
run_LIMA ${SAMPLES} ${TS_output}/CCS ${TS_output}/LIMA multiplex
run_REFINE ${SAMPLES} ${TS_output}/LIMA ${TS_output}/REFINE 
run_CLUSTER ${SAMPLES} ${TS_output}/REFINE ${TS_output}/CLUSTER
