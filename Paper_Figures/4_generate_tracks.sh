#!/bin/bash

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
OUTPUT_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/01_figures_tables/bash
rTg4510=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 

# subset bed file
mOntBed=${rTg4510}/G_Merged_Targeted/B_cupcake_pipeline/4_characterise/bed12Files/all_iso_ont_collapsed.filtered_counts_filtered_bothONTcounts_coloured.bed12
subset_fasta_gtf.py ${mOntBed} --bed --dir=${OUTPUT_DIR} --output=targetedNovelOntDiff -I=PB.14646.39352,PB.3948.4511,PB.14646.35283,PB.40586.1023
subset_fasta_gtf.py ${mOntBed} --bed --dir=${OUTPUT_DIR} --output=targetedProgOntDiff -I=PB.14646.39341,PB.14646.139,PB.20818.54,PB.20818.62,PB.14646.39352,PB.3948.4511,PB.14646.35283,PB.40586.1023,PB.41115.1365,PB.19309.7564,PB.40586.875,PB.20818.55