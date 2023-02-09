#!/bin/bash


##-------------------------------------------------------------------------
#************************************* Run TappAS on Knight
ISO_DIFF_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/3_differential/1_Input/A_IsoSeq_solo
HYBRID_DIFF_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole/3_differential/1_Input/B_hybrid_59

scp -r sl693@login.isca.ex.ac.uk:${ISO_DIFF_DIR}/* /mnt/data1/Szi/TAPPAS_FINAL_MOUSE/MouseWhole/IsoSeq_Expression
scp -r sl693@login.isca.ex.ac.uk:${ISO_DIFF_DIR}/* /mnt/data1/Szi/TAPPAS_FINAL_MOUSE/MouseWhole/RNASeq_Expression


##-------------------------------------------------------------------------
#************************************* Transfer results back to ISCA

### TappAS projects
# Project.030611158.tappas = Whole Mouse + IsoSeq
# Project.01501611467.tappas = Whole Mouse + RNASeq

SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome
source $SC_ROOT/2_Differential_Analysis/rTg4510_differential.config

cd ${TAPPAS_OUTPUT_DIR}; mkdir -p IsoSeq_Expression RNASeq_Expression

scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.030611158.tappas/* ${TAPPAS_OUTPUT_DIR}/IsoSeq_Expression/
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.01501611467.tappas/* ${TAPPAS_OUTPUT_DIR}/RNASeq_Expression/


##-------------------------------------------------------------------------
#************************************* Generate stats
Rscript $DIFFFUNC/IsoSeq_Whole_DEADIU.R $DiffAnalysis_WKD/TAPPAS_OUTPUT $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq_ISMrem.classification.txt $DiffAnalysis_WKD/TAPPAS_OUTPUT
cpat.py -x $REFERENCE/CPAT/Mouse_Hexamer.tsv -d $REFERENCE/CPAT/Mouse_logitModel.RData -g $DiffAnalysis_WKD/ORF/Cisd3.fasta --min-orf=50 --top-orf=50 -o Cisd3
