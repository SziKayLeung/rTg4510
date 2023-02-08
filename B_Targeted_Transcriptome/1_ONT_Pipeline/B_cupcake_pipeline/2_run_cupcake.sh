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
#SBATCH --output=2_run_cupcake.o
#SBATCH --error=2_run_cupcake.e

# 07/02/2023: Run minimap2, isoseq3 collapse and sqanti3 on ONT Targeted rTg4510 dataset

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/1_ONT_Pipeline/rTg4510_ont.config
source $SC_ROOT/B_Targeted_Transcriptome/1_ONT_Pipeline/01_source_functions.sh
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing 

# generate directories and specific variables
mkdir -p ${WKD_ROOT}/5_tclean/AllBatch2Batch3
mkdir -p ${WKD_ROOT}/6b_tofu_sqanti3
export MERGEID_DIR=${WKD_ROOT}/5_tclean/AllBatch2Batch3
export MERGEONT_DIR=${WKD_ROOT}/6b_tofu_sqanti3
export NAME=ONT_merged

# directory hierarchy
cd ${MERGEONT_DIR}
mkdir -p 1_minimap 2_collapse 3_sqanti3


##-------------------------------------------------------------------------
# Prepare for pipeline 
# 1. create AllBatch2Batch3_mergedsample_id.csv
# 2. create ONTBatch2Batch3_merged.fa

# generate _sample_id.csv for Batch 2 and Batch 3
source activate sqanti2_py3
adapt_cupcake_to_ont.py ${ONT_TCLEAN_DIR}/Batch2 -d ${MERGEID_DIR} -o AllBatch2 --batch
adapt_cupcake_to_ont.py ${ONT_TCLEAN_DIR}/Batch3 -d ${MERGEID_DIR} -o AllBatch3 --batch

# merge _sample_id.csv for Batch 2 and Batch 3
# head = keeps the header
# tail = take all the content from 2nd row downwards in sample_id.csv to merged file
head -n 1 ${MERGEID_DIR}/AllBatch3_sample_id.csv > ${MERGEID_DIR}/AllBatch2Batch3_mergedsample_id.csv
tail -n +2 -q ${MERGEID_DIR}/*_sample_id.csv >> ${MERGEID_DIR}/AllBatch2Batch3_mergedsample_id.csv

# merge fasta files from transcript clean for downstream alignment
print("Merging files:")
ls ${WKD_ROOT}/5_tclean/*/*.fa
tclean_fa=$(ls ${WKD_ROOT}/5_tclean/*/*.fa) 
cat ${tclean_fa} > ${MERGEID_DIR}/ONTBatch2Batch3_merged.fa


##-------------------------------------------------------------------------
# Run pipeline

source activate isoseq3
# align
cd ${MERGEONT_DIR}/1_minimap
pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${MERGEID_DIR}/ONTBatch2Batch3_merged.fa ${NAME}_mapped.bam \
  --unmapped --log-level TRACE --log-file ${NAME}_mapped.log

# collapse
cd ${MERGEONT_DIR}/2_collapse
isoseq3 collapse ${MERGEONT_DIR}/1_minimap/${NAME}_mapped.bam ${NAME}.gff --do-not-collapse-extra-5exons \
  --log-level TRACE --log-file ${NAME}_collapsed.log

# demultiplex
source activate nanopore
demux_cupcake_collapse.py ${MERGEONT_DIR}/2_collapse/${NAME}.read_stat.txt ${MERGEID_DIR}/AllBatch2Batch3_mergedsample_id.csv \
  --output=${NAME} --sample=${BARCODE_CONFIG} 

# sqanti3
cd ${MERGEONT_DIR}/3_sqanti3
python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf ${MERGEONT_DIR}/2_collapse/${NAME}.gff \
  $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report skip &> ${NAME}.sqanti.qc.log

# sqanti3 filter
python $SQANTI3_DIR/sqanti3_filter.py rules ${NAME}"_classification.txt" \
  --faa=${NAME}"_corrected.fasta" \
  --gtf=${NAME}"_corrected.gtf" \
  -j=${SQANTIFIL_JSON} \
  --skip_report &> ${NAME}.sqanti.filter.log
