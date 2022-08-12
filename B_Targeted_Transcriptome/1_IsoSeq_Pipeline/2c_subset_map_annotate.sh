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
#SBATCH --output=Targeted_Mouse_Part3.o
#SBATCH --error=Targeted_Mouse_Part3.e

# 03/06/2021: Same script as Targeted_Mouse_Part2.sh but only work with the same samples as Whole Transcriptome for fair comparisons
# 14/01/2022: Re-run SQANIT3 on subset of whole samples with no RNA-Seq Support
##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh

##-------------------------------------------------------------------------
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch: Functions 1 - 3, 6 for all batched runs"
echo "Already processed in individual samples for full set: Functions 4"

##-------------------------------------------------------------------------
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
# merging_at_refine <output_name> <output_root_dir> <samples.....>
merging_at_refine $SUBNAME $WKD_ROOT $SUBSET_DIR ${SUBSET_SAMPLES_NAMES[@]}

# run_map_cupcakecollapse <output_name> <root_input_directory> <root_output_directory>
run_map_cupcakecollapse $SUBNAME $WKD_ROOT $SUBSET_DIR

# demux_targeted <output_name> <input_root_dir> <output_root_dir>
demux_targeted $SUBNAME $WKD_ROOT $SUBSET_DIR

##-------------------------------------------------------------------------
echo "#************************************* RNAseq [Function 9, 10]"
echo "Already processed in individual samples for rnaseq alignment for full set: Functions 9,10"

##-------------------------------------------------------------------------
echo "#************************************* RNASeq & IsoSeq [Function 11]"
# run_kallisto <sample_prefix_output_name> <io_dir> 
run_kallisto $SUBNAME $SUBSET_DIR 

##-------------------------------------------------------------------------
echo "#************************************* SQANTI3 [Function 12]"
## 12) run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna> <io_dir>
run_sqanti3 $SUBNAME full $SUBSET_DIR
run_sqanti3 $SUBNAME basic $SUBSET_DIR

##-------------------------------------------------------------------------
echo "#************************************* TAMA filter [Function 13,14,16]"
## 13) TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_remove_fragments $SUBNAME full

## 14) TAMA_sqanti_filter <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_sqanti_filter $SUBNAME full

## 16) TAMA_tappas_input <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_tappas_input $SUBNAME full


