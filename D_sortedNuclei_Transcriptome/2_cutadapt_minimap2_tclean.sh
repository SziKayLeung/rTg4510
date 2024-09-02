#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-7%4 # 8 barcodes per flow cell
#SBATCH --output=Output/2log/2_cutadapt_minimap2_tclean-%A_%a.o
#SBATCH --error=Output/2log/2_cutadapt_minimap2_tclean-%A_%a.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/D_sortedNuclei_Transcriptome
source $SC_ROOT/rTg4510_snont.config
source $SC_ROOT/01_source_functions.sh

NeuNsample=${NEUN_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
DNsample=${DN_SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}

##-------------------------------------------------------------------------

# merge each sample into one fastq file 
merge_fastq_across_samples ${NeuNsample} ${WKD_ROOT}/1_demultiplex/NeuN ${WKD_ROOT}/1b_demultiplex_merged/NeuN
merge_fastq_across_samples ${DNsample} ${WKD_ROOT}/1_demultiplex/DN ${WKD_ROOT}/1b_demultiplex_merged/DN

# delinate polyA and polyT sequences, reverse complement polyT sequences, remove polyA from all sequences
post_porechop_run_cutadapt ${WKD_ROOT}/1b_demultiplex_merged/NeuN/${NeuNsample}_merged.fastq ${WKD_ROOT}/2_cutadapt_merge/NeuN
post_porechop_run_cutadapt ${WKD_ROOT}/1b_demultiplex_merged/DN/${DNsample}_merged.fastq ${WKD_ROOT}/2_cutadapt_merge/DN
 
# map combined fasta to reference genome
run_minimap2 ${WKD_ROOT}/2_cutadapt_merge/NeuN/${NeuNsample}_merged_combined.fasta ${WKD_ROOT}/3_minimap/NeuN
run_minimap2 ${WKD_ROOT}/2_cutadapt_merge/DN/${DNsample}_merged_combined.fasta ${WKD_ROOT}/3_minimap/DN

# run transcript clean on aligned reads
run_transcriptclean ${WKD_ROOT}/3_minimap/NeuN/${NeuNsample}_merged_combined_sorted.sam ${WKD_ROOT}/4_tclean/NeuN
run_transcriptclean ${WKD_ROOT}/3_minimap/DN/${DNsample}_merged_combined_sorted.sam ${WKD_ROOT}/4_tclean/DN