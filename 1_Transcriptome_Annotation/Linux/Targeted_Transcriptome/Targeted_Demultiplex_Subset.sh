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
#SBATCH --output=Targeted_Demutliplex_Subset.o
#SBATCH --error=Targeted_Demutliplex_Subset.e

# 20/01/2021: IsoSeq Analysis pipeline for merging individual samples that were also sequenced as whole for comparison
  # No need for demultiplexing, only proceed from merge at refine cluster onwards
# 21/01/2021: Alternative pipeline for targeted

#************************************* DEFINE GLOBAL VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
DEMUX_SCRIPT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome/
SAMPLES_LIST=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/cupcake/tofu/counting
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq/REFINE
All_Samples=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged_Subset
Targeted_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/
cd $All_Samples; mkdir CLUSTER MAPPING TOFU SQANTI

#************************************* All samples, WT, TG merged at refine
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
source $FUNCTIONS/Post_Isoseq3_Function.sh
module load Miniconda2
module load R

#run_pipeline_to_chain <Sample_Name> <List.of.Samples...>
run_pipeline_to_tofu(){

  SAMPLES=$(echo "${@:2}")
  echo "Processing together: $SAMPLES"

  # merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
  # convert_fa2fq <file_name> <input_dir>
  # run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
  # tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir>
  # demux <input path read.stat file> <input path of samples file> <path of output>
    # run Cupcake_Demultiplex.R, read in read.stat file from cupcake collapse output and count abundance of each sample based on CCS_ID

  merging_at_refine $REFINE $All_Samples/CLUSTER $1 ${SAMPLES[@]}
  convert_fa2fq $1".clustered.hq.fasta" $All_Samples/CLUSTER
  run_minimap2 $1 $All_Samples/CLUSTER mm10 $All_Samples/MAPPING
  on_target_rate $Targeted_dir/Probes/FINAL_MOUSE.bed $All_Samples/CLUSTER/All_Targeted_Merged_Subset.clustered.hq.fasta $All_Samples/MAPPING/All_Targeted_Merged_Subset.clustered.hq.fastq.sam $All_Samples/MAPPING/All_Targeted_Merged_Subset.fasta.sam.probe_hit.txt
  tofu $1 $All_Samples/CLUSTER $All_Samples/MAPPING $All_Samples/TOFU

  # demultiplex collapsed.read_stat.txt
  Rscript $DEMUX_SCRIPT/Demultiplex_Cupcake.R
}

run_pipeline_to_tofu All_Targeted_Merged_Subset O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23
# run_sqanti2_QC <sample_name> <working_directory> <collapsed.gtf> <count.txt>
run_sqanti2 All_Targeted_Merged_Subset $All_Samples/SQANTI $All_Samples/TOFU/All_Targeted_Merged_Subset.collapsed.filtered.gff $All_Samples/TOFU/All_Targeted_Merged_Subset.Demultipled_Abundance.txt
