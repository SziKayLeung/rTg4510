#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=30:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-1 # 2 samples
#SBATCH --output=TargetedONT-%A_%a.o
#SBATCH --error=TargetedONT-%A_%a.e

### SBATCH -p mrchq  #SBATCH --mem=500G ### for porechop

# Aim: Run Batch 2 and Batch 3 simultaneoulsy: nanofilt, porechop, post-targeted-porechop, cutadapt and combine
# Dataset: ONT Targeted rTg4510 dataset (n = 18 samples)

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/1_ONT_Pipeline/rTg4510_ont.config
source $SC_ROOT/B_Targeted_Transcriptome/1_ONT_Pipeline/01_source_functions.sh

##-------------------------------------------------------------------------

# run as array (defined in config file)
BATCH=${BATCH_NAMES[${SLURM_ARRAY_TASK_ID}]}
FASTA=${FASTA_FILES[${SLURM_ARRAY_TASK_ID}]}


##-------------------------------------------------------------------------

echo "Running pipeline for ${BATCH}"
# 1) run_merge <input_fasta_dir> <sample_output_name>
run_merge ${FASTA} ${BATCH}

# 2) run_nanofilt <sample_output_name_merged.fastq> <sample_output_name> <output_directory>
run_nanofilt $RAW_ROOT_DIR/${BATCH}"_Merged.fq" ${BATCH} $WKD_ROOT

# 3) run_porechop <sample> <root_directory> <type>
run_porechop ${BATCH} $WKD_ROOT Targeted

# 4) post_porechop <sample> <root_directory>
post_targeted_porechop ${BATCH} $WKD_ROOT

# 5) run_cutadapt_and_combine <sample> <root_directory> <demux/nondemux>
run_cutadapt_and_combine ${BATCH} $WKD_ROOT nondemux

# 6) run_minimap2 <sample> <root_directory> <batch>
# 7) run_transcriptclean <sample> <root_dir> <batch>
# 8) run_talon_label <sample> <root_dir> <batch>
for demuxsample in ${ALL_BARCODE_NAMES[@]}; do
  
  run_minimap2 $demuxsample $WKD_ROOT ${BATCH}
  run_transcriptclean $demuxsample $WKD_ROOT ${BATCH}
  run_talon_label $demuxsample $WKD_ROOT ${BATCH}
  
done
  