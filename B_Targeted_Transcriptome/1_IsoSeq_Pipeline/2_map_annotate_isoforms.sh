#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Targeted_ADBDR_Part2c.o
#SBATCH --error=Targeted_ADBDR_Part2c.e

# 11/05/2021: Run pipeline from refine to sqanti and tama (output: Targeted_ADBDR_Part2.o)
# 03/11/2021: Run kallisto with ADBDR RNA-Seq (30 samples), SQANTI3 with RNA-Seq support and TAMA filter (output: Targeted_ADBDR_Part2b.o)

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch (1_run_isoseq3.sh) for all batched runs"

## 4) run_targeted_REFINE 
run_targeted_REFINE 

# 6) run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
for i in ${ALL_SAMPLES_NAMES[@]}; do run_CLUSTER $i; done

##-------------------------------------------------------------------------
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
# 5) merging_at_refine <output_name> <output_root_dir> <samples.....>
merging_at_refine $NAME $WKD_ROOT $WKD_ROOT ${ALL_SAMPLES_NAMES[@]}

# 7) run_map_cupcakecollapse <output_name> <root_input_directory> <root_output_directory>
run_map_cupcakecollapse $NAME $WKD_ROOT $WKD_ROOT

# demux_targeted <output_name> <input_root_dir> <output_root_dir>
demux_targeted $SUBNAME $WKD_ROOT $WKD_ROOT

##-------------------------------------------------------------------------
echo "#************************************* RNAseq [Function 9, 10]"
## 9) run_star <list_of_samples> <input_directory> <output_dir>
bash $SC_ROOT/1_IsoSeq_Pipeline/2b_map_rnaseq_genome.sh; wait

## 10) mouse_merge_fastq <sample>
mouse_merge_fastq $NAME

##-------------------------------------------------------------------------
echo "#************************************* RNASeq & IsoSeq [Function 11]"
# run_kallisto <sample_prefix_output_name> <io_dir> 
run_kallisto $NAME $WKD_ROOT

##-------------------------------------------------------------------------
echo "#************************************* SQANTI3 [Function 12]"
## 12) run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna/nornaseq>
run_sqanti3 $NAME full
run_sqanti3 $NAME basic
run_sqanti3 $NAME nornaseq

##-------------------------------------------------------------------------
echo "#************************************* TAMA filter [Function 13,14,16]"
## 13) TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna/nornaseq> <io_dir>
TAMA_remove_fragments $NAME full $WKD_ROOT

## 14) TAMA_sqanti_filter <sample> <mode=basic/full/nokallisto/lncrna> <io_dir>
TAMA_sqanti_filter $NAME full $WKD_ROOT

## 16) TAMA_tappas_input <sample> <mode=basic/full/nokallisto/lncrna> <io_dir>
TAMA_tappas_input $NAME full $WKD_ROOT

# remove_3ISM <sample> <mode=basic/full/nokallisto/lncrna/nornaseq>
remove_3ISM $NAME nornaseq

##-------------------------------------------------------------------------
echo "#************************************* QC [Function 15]"
# parse_stats_per_sample <sample>
parse_stats_per_sample $NAME 

run_target_rate

##-------------------------------------------------------------------------
echo "#************************************* Characterisation of isoforms"
find_humanMAPT

##-------------------------------------------------------------------------
echo "#************************************* Whole vs Targeted Transcriptome "
# In Targeted_Part3.sh