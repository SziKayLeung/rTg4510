#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=15:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --reservation=research_project-mrc148213_5

#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis_noRNASEQ
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/RNASeq

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome
GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
TAPPASFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis

module load Miniconda2/4.3.21
dataset=AllMouseTargeted

##################################################################################################
#************************************* Prepare files for TappAS
# 7 files:
# 1) Iso-Seq Annotation scaffold file                       = WholeIsoSeq.collapsed.gff3
# 2) Iso-Seq retained isoforms from tama filtering          = tama_retained_pbid.txt
# 3) Iso-Seq Expression File after further collapse
# 4) Iso-Seq Phenotype File
# 5) RNA-Seq Expression File after alignment to Iso-Seq collapsed by length
# 6) RNA-Seq Phenotype File


# Targeted Transcriptome: RNA-Seq alignment to Whole + Targeted Transcriptome
# 8) RNA-Seq Expression file - Kallisto RNA-Seq alignment to tama merged transcriptome of whole and targeted (all isoforms)
# 9) Annotation file - SQANTI3 annotation from tama merged transcriptome
# 10) List of Isoforms from target genes

# File 1,2
cp $DiffAnalysis_WKD/SQANTI3/$dataset".collapsed.gff3" $DiffAnalysis_WKD/TAPPAS_INPUT
cp $DiffAnalysis_WKD/SQANTI_TAMA_FILTER/tama_retained_pbid.txt $DiffAnalysis_WKD/TAPPAS_INPUT

# File 3, 4, - using Iso-Seq Expression
# counts_subset_4tappas <input_class> <output_class> <type_genes>
# regenerate expression matrix for tappas from long reads FL read counts
counts_subset_4tappas $DiffAnalysis_WKD/COLLAPSE_FILTER/$dataset"_sqantisubset.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression/$dataset"_sqantisubset.expression.txt" AD
cp $RAWDIR/TargetedMouse_PhenotypeTAPPAS.txt $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression

# File 5, 6 - using RNA-Seq Expression
# Rscript script.R <input.dir> <output.file> <type=Whole/Targeted/WholeTargeted> <targeted.class.files>
source activate sqanti2_py3; Rscript $TAPPASFUNC/TAPPAS_RNASEQ_Exp.R $DiffAnalysis_WKD/RNASeq_SQANTI3 TargetedMouseRNASeq.expression.txt Targeted $DiffAnalysis_WKD/COLLAPSE_FILTER/$dataset"_sqantisubset.classification.txt"
cp $DiffAnalysis_WKD/RNASeq_SQANTI3/TargetedMouseRNASeq.expression.txt $DiffAnalysis_WKD/TAPPAS_INPUT/RNASeq_Expression
cp $RAWDIR/TargetedMouse_RNASeqPhenotypeTAPPAS.txt $DiffAnalysis_WKD/TAPPAS_INPUT/RNASeq_Expression


##################################################################################################
#************************************* Run TappAS on Knight
scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/TAPPAS_INPUT/* /mnt/data1/Szi/TAPPAS_MouseTargetedSQ3
/mnt/data1/Aaron/sw/jre1.8.0_181/bin/java -jar /mnt/data1/Szi/tappAS.jar
