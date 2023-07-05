#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=50:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --error=3_merged_characterise.e
#SBATCH --output=3_merged_characterise.o

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
## Characterisation with CPAT and colour by abundance

# LOGEN: subset cupcake classification file by target genes
# merge cupcake classification file with abundance
source activate nanopore
subset_quantify_filter_tgenes.R \
--classfile ${CUPMERGE_DIR}/3_sqanti3/${MERGED_NAME}_collapsed_RulesFilter_result_classification.txt \
--expression ${CUPMERGE_DIR}/2_collapse/demux_fl_count.csv \
--target_genes ${TGENES_TXT} 

# filter cupcake classification file with minimum number of reads and counts
subset_quantify_filter_tgenes.R \
--classfile ${CUPMERGE_DIR}/3_sqanti3/${MERGED_NAME}_collapsed_RulesFilter_result_classification.txt \
--expression ${CUPMERGE_DIR}/2_collapse/demux_fl_count.csv \
--target_genes ${TGENES_TXT} \
--filter --nsample=5 --nreads=10


# working variables
finalanno=${CUPMERGE_DIR}/3_sqanti3/${MERGED_NAME}_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt 
finaliso=${CUPMERGE_DIR}/3_sqanti3/${MERGED_NAME}_collapsed_RulesFilter_result_classification.targetgenes_filtered_isoforms.txt

# LOGEN: subset fasta and gtf using the finalised list of target gene isoforms
source activate sqanti2_py3
subset_fasta_gtf.py --gtf ${CUPMERGE_DIR}/3_sqanti3/all_iso_ont_collapsed.filtered.gtf -i ${finaliso} -o counts_filtered
subset_fasta_gtf.py --fa ${CUPMERGE_DIR}/3_sqanti3/all_iso_ont_collapsed.filtered.faa -i ${finaliso} -o counts_filtered

# run_cpat <input_fasta> <output_name> <output_dir>
run_cpat ${CUPMERGE_DIR}/3_sqanti3/${MERGED_NAME}_collapsed.filtered_counts_filtered.fa $MERGED_NAME ${CUPMERGE_DIR}

# extract_best_orf <sample> <root_dir>
extract_best_orf $MERGED_NAME ${CUPMERGE_DIR}

# colour_by_abundance <sample> <input_gtf> <abundance_file> <root_dir>
colour_by_abundance ${MERGED_NAME} \
  ${CUPMERGE_DIR}/3_sqanti3/${MERGED_NAME}"_collapsed.filtered_counts_filtered.gtf" \
  ${CUPMERGE_DIR}/2_collapse/demux_fl_count.csv ${CUPMERGE_DIR}

# run_transdecoder <name> <root_dir>
run_transdecoder ${MERGED_NAME} ${CUPMERGE_DIR}

##-------------------------------------------------------------------------
## Characterisation with FICLE

# FICLE: subset_gene_reference 
subset_gene_reference ${CUPMERGE_DIR}

# FICLE: full characterisation
# coloured by ONT and Iso-Seq
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

# FICLE: extract reference information
extract_reference_info.py --r ${GENOME_GTF} --glist ${TGENES[@]} \
  --split ${CUPMERGE_DIR}/4_characterise/TargetGenesRef \
  --short_read  ${RNASEQ_COUNTS} 