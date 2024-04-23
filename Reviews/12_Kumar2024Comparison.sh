#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=12_Kumar2024Comparison.o
#SBATCH --error=12_Kumar2024Comparison.o

# 23/04/2024: Run Kumar dataset through ONT pipeline for downstream comparison

##-------------------------------------------------------------------------

# https://academic.oup.com/nar/article/52/6/2865/7627472#supplementary-data
# https://figshare.com/articles/dataset/Long_read_oxford_nanopore_technology_mouse_brain_aging_datasets/25670262

export K2024=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/0_otherStudies/Kumar2024
export REFERENCE=/lustre/projects/Research_Project-MRC148213/lsl693/references
export GENOME_FASTA=$REFERENCE/mouse/mm10.fa
export GENOME_GTF=$REFERENCE/annotation/gencode.vM22.annotation.gtf
export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
export CUPCAKE=$SOFTDIR/cDNA_Cupcake
export ANNOTATION=$CUPCAKE/annotation
export SEQUENCE=$CUPCAKE/sequence
export PYTHONPATH=$PYTHONPATH:$SEQUENCE
export SQANTI3_DIR=$SOFTDIR/SQANTI3
export SQANTIFIL_JSON=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/0_utils/filter_default_reducecoverage.json
export NAME=Kumar2024

module load Miniconda2/4.3.21
source activate isoseq3


##-------------------------------------------------------------------------
# already aligned, convert to fasta
#cd ${K2024}/aligned
#for i in ${K2024}/aligned/*bam*; do 
#  echo $i
#  sample=$(basename "$i" | cut -d "." -f 1 )
#  echo $sample
#  samtools bam2fq $i| seqtk seq -A > $sample.fa
#done

#aligned=$(ls ${K2024}/aligned/*.fa) 
#cat ${aligned} > merged.fa


##-------------------------------------------------------------------------

# isoseq collapse
cd ${K2024}/cupcake
#pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${K2024}/aligned/merged.fa ${NAME}_mapped.bam --unmapped --log-level TRACE --log-file ${NAME}_mapped.log
#isoseq3 collapse ${NAME}_mapped.bam ${NAME}.gff --do-not-collapse-extra-5exons --log-level TRACE --log-file ${NAME}"_collapsed.log"

# isoseq sqanti3
source activate sqanti2_py3
python $SQANTI3_DIR/sqanti3_qc.py -t 30 ${NAME}.gff $GENOME_GTF $GENOME_FASTA --genename --skipORF --report skip &> merged.sqanti.qc.log
python $SQANTI3_DIR/sqanti3_filter.py rules ${NAME}"_classification.txt" --faa=${NAME}"_corrected.fasta" --gtf=${NAME}"_corrected.gtf" \
  -j=${SQANTIFIL_JSON} --skip_report &> ${NAME}.sqanti.filter.log
  
rTg4510_gtf=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/B_cupcake_pipeline/3_sqanti3/all_iso_ont_collapsed.filtered_counts_filtered.gtf
cd ${K2024}/comparison
cp ${rTg4510_gtf} .
cp ${K2024}/cupcake/Kumar2024.filtered.gtf .
PATH="/lustre/projects/Research_Project-MRC148213/lsl693/software/gffcompare:$PATH"
gffcompare -r all_iso_ont_collapsed.filtered_counts_filtered.gtf Kumar2024.filtered.gtf -o rTg4510Kumar
gffcompare -r Kumar2024.filtered.gtf all_iso_ont_collapsed.filtered_counts_filtered.gtf -o KumarrTg4510
