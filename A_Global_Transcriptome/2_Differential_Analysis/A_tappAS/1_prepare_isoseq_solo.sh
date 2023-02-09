#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

## ---------------------------
## Purpose: Generate files (Iso-Seq annotation and expression) needed for tappAS 
## 
## ---------------------------

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/2_Differential_Analysis/01_source_functions.sh
source $SC_ROOT/2_Differential_Analysis/rTg4510_differential.config


##-------------------------------------------------------------------------
#************************************* Prepare files for TappAS
# 4 files:
# 1) Iso-Seq Annotation scaffold file                       = XX.collapsed.gff3
# 2) Iso-Seq retained isoforms                              
# 3) Iso-Seq Expression File
# 4) Iso-Seq Phenotype File


# create directory 
mkdir -p ${TAPPAS_INPUT_DIR}/A_IsoSeq_solo

# File 1
# Annotation file generated from IsoAnnotLite in SQANTI3
# Note, RNA-Seq reads were not used as junction filter, given lower depth and reduce false negative
cp $ISO_WKD_ROOT/${NAME}.collapsed.gff3 ${TAPPAS_INPUT_DIR}/A_IsoSeq_solo

# File 2
cp $ISO_WKD_ROOT/${NAME}_ISMrem.isoform.txt ${TAPPAS_INPUT_DIR}/A_IsoSeq_solo

# File 3
source activate nanopore
Rscript ${PREPCOUNTS} $ISO_WKD_ROOT/${NAME}.collapsed_classification.filtered_lite_classification.txt ${TAPPAS_INPUT_DIR}/A_IsoSeq_solo/${NAME}_expression.txt

# File 4
cp $ISO_TAPPAS_PHENO ${TAPPAS_INPUT_DIR}/A_IsoSeq_solo


