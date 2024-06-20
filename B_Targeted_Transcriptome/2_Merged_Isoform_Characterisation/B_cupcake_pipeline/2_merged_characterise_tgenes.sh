#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --error=2_merged_characterise_tgenes.e
#SBATCH --output=2_merged_characterise_tgenes.o

# 13/02/2023: characterisation of merged Iso-Seq + ONT transcripts with FICLE

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/rTg4510_merged.config
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/01_source_function.sh
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
export PATH=$PATH:${LOGEN_ROOT}/merge_characterise_dataset
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 
export PATH=$PATH:${FICLE_ROOT}
export PATH=$PATH:${FICLE_ROOT}/reference


##-------------------------------------------------------------------------

# FICLE: full characterisation
# coloured by ONT and Iso-Seq
source activate sqanti2_py3
for g in ${TGENES[@]}; do

  echo $g
  mkdir -p ${CUPMERGE_DIR}/4_characterise/TargetGenes
  mkdir -p ${CUPMERGE_DIR}/4_characterise/TargetGenes/Log
  output_dir=${CUPMERGE_DIR}/4_characterise/TargetGenes
  
  ficle.py -n=$g \
  -r=${CUPMERGE_DIR}/4_characterise/TargetGenesRef/ \
  -b=${CUPMERGE_DIR}/4_characterise/bed12Files/${MERGED_NAME}_collapsed.filtered_counts_filtered_concat_counts_coloured.bed12 \
  -g=${CUPMERGE_DIR}/3_sqanti3/${MERGED_NAME}_collapsed.filtered_counts_filtered.gtf \
  -c=${CUPMERGE_DIR}/3_sqanti3/${MERGED_NAME}_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt \
  --cpat=${CUPMERGE_DIR}/4_characterise/CPAT/${MERGED_NAME}.ORF_prob.best.tsv   \
  -o=$output_dir &> ${CUPMERGE_DIR}/4_characterise/TargetGenes/Log/$g"_characterise.log"

done