#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=../../../bash_output/1c_merged_collapse_sqanti3.o
#SBATCH --error=../../../bash_output/1c_merged_collapse_sqanti3.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/rTg4510_merged.config
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/01_source_function.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous

export dir=${CUPMERGE_DIR}
export samplename=all_iso_ont

##-------------------------------------------------------------------------

# merge alignment
echo "Collapsing..."
allfilteredmapped=($(ls ${dir}/1_align/*mapped.filtered.sorted.bam)) 
ls ${allfilteredmapped[@]}
source activate nanopore
samtools merge -f ${dir}/2_collapse/${samplename}_mapped.filtered.sorted.bam ${allfilteredmapped[@]}

# collapse
echo "Collapsing..."
echo "Output: ${dir}/2_collapse/${samplename}_collapsed.gff"
cd ${dir}/2_collapse
source activate isoseq3
isoseq3 collapse ${dir}/2_collapse/${samplename}_mapped.filtered.sorted.bam ${samplename}_collapsed.gff \
  --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
  --log-level TRACE --log-file ${samplename}_collapsed.log

# demultiplex 
# demux_ont_isoseq_cupcake_collapse.py <merged_read_stat.txt> <ont_sample_id.csv> <iso_sample_id.csv>
source activate sqanti2_py3
demux_ont_isoseq_cupcake_collapse.py \
  ${dir}/2_collapse/${samplename}_collapsed.read_stat.txt \
  ${ONT_TCLEAN_DIR}/AllBatch2Batch3/AllBatch2Batch3_mergedsample_id.csv \
  ${ISO_MERGED_CLUSTER_DIR}/AllMouseTargeted.clustered.demuxed.cluster_report.csv \
  --ont_sample=${ONT_BARCODE_CONFIG} 

# sqanti3
echo "Running SQANTI3..."
cd ${dir}/3_sqanti3
  python $SQANTI3_DIR/sqanti3_qc.py ${dir}/2_collapse/${samplename}_collapsed.gff \
  $GENOME_GTF $GENOME_FASTA -t 30 --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report skip &> ${samplename}_sqanti_qc.log

# sqanti3 filter 
filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_collapsed_classification.txt \
  --faa=${samplename}_collapsed_corrected.fasta \
  --gtf=${samplename}_collapsed_corrected.gtf \
  -j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log