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


##-------------------------------------------------------------------------

# input 
module load Miniconda2
LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/D_sortedNuclei_Transcriptome
source $SC_ROOT/rTg4510_snont.config
source $SC_ROOT/01_source_functions.sh

sortedBam=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/5_cupcake/6_collapse/rTg4510SCN_mapped.filtered.sorted.bam
bulkBam=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/B_cupcake_pipeline/2_collapse/all_iso_ont_mapped.filtered.sorted.bam
sortedDir=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/5_cupcake/7_sqanti3
bulkDir=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/B_cupcake_pipeline/3_sqanti3


##-------------------------------------------------------------------------

# merge bam files
#source activate nanopore
#cd /lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged
#samtools merge -f bulkSorted_mapped.filtered.sorted.bam ${sortedBam} ${bulkBam}

# run isoseq3
#source activate isoseq3
#isoseq3 collapse bulkSorted_mapped.filtered.sorted.bam bulkSorted_collapsed.gff \
#--min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
#--log-level TRACE --log-file bulkSorted_collapsed.log


## ----- run gffcompare ----

#source activate sqanti2_py3

#cd /lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged

# extract gtf from target genes in sorted nuclei data
#awk '{print $1}' ${sortedDir}/rTg4510SCN_collapsed_RulesFilter_result_classification_targetgenes_counts.txt > ${sortedDir}/rTg4510SCN_collapsed_RulesFilter_result_classification_targetgenes_isoforms.txt
#subset_fasta_gtf.py --gtf ${sortedDir}/rTg4510SCN_collapsed.filtered.gtf -i ${sortedDir}/rTg4510SCN_collapsed_RulesFilter_result_classification_targetgenes_isoforms.txt -o targetgenes_filtered


#sortedGtf=${sortedDir}/rTg4510SCN_collapsed.filtered_targetgenes_filtered.gtf
#bulkGtf=${bulkDir}/all_iso_ont_collapsed.filtered_counts_filtered.gtf

#cp ${sortedGtf} .
#cp ${bulkGtf} .
#PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/gffcompare:$PATH"
#gffcompare -r ${bulkGtf} ${sortedGtf} -o bulkSorted
#gffcompare -r ${sortedGtf} ${bulkGtf} -o Sortedbulk

#cp ${sortedDir}/*map* .
#cp ${bulkDir}/*map* .


##-------------------------------------------------------------------------

# merge fasta files
source activate nanopore
outputDir=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged
cd /lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/6_merged
#subset_fasta_gtf.py --fa ${sortedDir}/rTg4510SCN_collapsed_corrected.fasta -i ${sortedDir}/rTg4510SCN_collapsed_RulesFilter_result_classification_targetgenes_isoforms.txt -o targetgenes_filtered

sortedFa=${sortedDir}/rTg4510SCN_collapsed_corrected_targetgenes_filtered.fa
bulkFa=${bulkDir}/all_iso_ont_collapsed.filtered_counts_filtered.fa
cat ${sortedFa} ${bulkFa} > merged.fa

source activate isoseq3
pbmm2 align --preset ISOSEQ --sort $GENOME_FASTA merged.fa merged_mapped.bam --log-level TRACE --log-file merged.log

filter_alignment merged_mapped ${outputDir}

source activate isoseq3
isoseq3 collapse merged_mapped.filtered.bam merged_collapse.gff --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons --log-level TRACE --log-file merged_collapsed.log

mkdir original_fasta
cp ${sortedFa} .
cp  ${bulkFa} .
adapt_cupcake_to_ont.py ${outputDir}/original_fasta -o merged 

source activate sqanti2_py3
python $SQANTI3_DIR/sqanti3_qc.py merged_collapse.gff \
$GENOME_GTF $GENOME_FASTA -t 30 --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA \
--genename --isoAnnotLite --skipORF --report skip


python $SQANTI3_DIR/sqanti3_filter.py rules merged_collapse_classification.txt \
--faa=merged_collapse_corrected.fasta \
--gtf=merged_collapse_corrected.gtf \
-j=${filteringJson} --skip_report &> merged_collapse_sqanti_filter.log

# run isoseq3
#source activate isoseq3
#isoseq3 collapse bulkSorted_mapped.filtered.sorted.bam bulkSorted_collapsed.gff \
#--min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
#--log-level TRACE --log-file bulkSorted_collapsed.log
