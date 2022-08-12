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
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq

cd $DiffAnalysis_WKD; mkdir -p TAPPAS_OUTPUT
cd $DiffAnalysis_WKD/TAPPAS_OUTPUT; mkdir -p IsoSeq_Expression RNASeq_Expression

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome
GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
TAPPASFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Rscripts/
DIFFFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Figures_Thesis/DiffAnalysis_WholeTranscriptome

module load Miniconda2/4.3.21
source activate sqanti2_py3

dataset=WholeIsoSeq

##################################################################################################
#************************************* Prepare files for TappAS
# 5 files:
# 1) Iso-Seq Annotation scaffold file                       = WholeIsoSeq.collapsed.gff3
# 2) Iso-Seq Retained file
# 3) Iso-Seq Expression File
# 4) Iso-Seq Phenotype File
# 5) RNA-Seq Expression File
# 6) RNA-Seq Phenotype File


# File 1
cp $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq.collapsed.gff3 $DiffAnalysis_WKD/TAPPAS_INPUT

# File 2
cp $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq_ISMrem.isoform.txt $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression

# File 3
Rscript $GENERALFUNC/Whole_Prepare_Counts.R $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression/WholeIsoSeq_expression.txt

# File 4
cp $RAWDIR/WholeIsoSeq_PhenotypeTAPPAS.txt $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression

# File 5
# Rscript script.R <input.dir> <output.file> <type=Whole/Targeted> <targeted.class.files>
Rscript $TAPPASFUNC/TAPPAS_RNASEQ_Exp.R $DiffAnalysis_WKD/RNASeq_SQANTI3 WholeMouseRNASeq_sqantisubset.expression.txt Whole NA
cp $DiffAnalysis_WKD/RNASeq_SQANTI3/WholeMouseRNASeq_sqantisubset.expression.txt $DiffAnalysis_WKD/TAPPAS_INPUT/RNASeq_Expression

# File 6
cp $RAWDIR/WholeAllMouse_PhenotypeTAPPAS.txt $DiffAnalysis_WKD/TAPPAS_INPUT/RNASeq_Expression

##################################################################################################
#************************************* Run TappAS on Knight
scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/TAPPAS_INPUT/* /mnt/data1/Szi/TAPPAS_FINAL_MOUSE/MouseWhole
/mnt/data1/Aaron/sw/jre1.8.0_181/bin/java -jar /mnt/data1/Szi/tappAS.jar

scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final/TAPPAS_INPUT/IsoSeq_Expression/WholeIsoSeq_expression.txt /mnt/data1/Szi/TAPPAS_FINAL_MOUSE/MouseWhole

#### Tappas output
### tappAS projects
# Project.030611158.tappas = Whole Mouse + IsoSeq
# Project.01501611467.tappas = Whole Mouse + RNASeq

scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.01501611467.tappas/* $DiffAnalysis_WKD/TAPPAS_OUTPUT/RNASeq_Expression/
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.030611158.tappas/* $DiffAnalysis_WKD/TAPPAS_OUTPUT/IsoSeq_Expression/

##################################################################################################
#************************************* Generate stats
Rscript $DIFFFUNC/IsoSeq_Whole_DEADIU.R $DiffAnalysis_WKD/TAPPAS_OUTPUT $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq_ISMrem.classification.txt $DiffAnalysis_WKD/TAPPAS_OUTPUT

# Others
grep PB.2833.1 -A 20 $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_ISMrem.fasta
grep PB.2833.2 -A 20 $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_ISMrem.fasta
grep PB.16934.2 -A 20 $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_ISMrem.fasta

grep PB.11626.2  -A 50 $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_ISMrem.fasta
grep PB.8363.3  -A 50 $DiffAnalysis_WKD/SQANTI3/WholeIsoSeq.collapsed_classification.filtered_lite_ISMrem.fasta


cpat.py -x $REFERENCE/CPAT/Mouse_Hexamer.tsv -d $REFERENCE/CPAT/Mouse_logitModel.RData -g $DiffAnalysis_WKD/ORF/Cisd3.fasta --min-orf=50 --top-orf=50 -o Cisd3
