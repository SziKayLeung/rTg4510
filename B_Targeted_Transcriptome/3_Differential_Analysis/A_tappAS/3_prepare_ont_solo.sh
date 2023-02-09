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
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON
SQANTI3_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON/MissingBC1/SQANTI3_Unfiltered_RNASeq
cd $DiffAnalysis_WKD; mkdir -p TAPPAS_INPUT TAPPAS_OUTPUT RNASeq_SQANTI3

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Nanopore/Targeted_Transcriptome
GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
TAPPASFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis
NANOPOREFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Nanopore

module load Miniconda2/4.3.21
source activate nanopore
dataset=AllMouseTargeted
sqname=$dataset".collapsed_classification.filtered_lite"
mergedsqpath=$DiffAnalysis_WKD/All/Merged/SQANTI3/IsoSeqONT_final_genename

##################################################################################################
#************************************* Prepare files for TappAS
# 4 files:
# 1) Iso-Seq retained isoforms from tama filtering          = tama_retained_pbid.txt
# 2) Iso-Seq Expression File after further collapse
# 3) Iso-Seq Annotation scaffold file                       = WholeIsoSeq.collapsed.gff3
# 4) Iso-Seq Phenotype File


# File 1, 2
# Note: ONT_unfiltered_expression.txt generated from XXX
# Note ONT_retained_list is generated from the Rscript Nanopore variables
Rscript $NANOPOREFUNC/Rscripts/normalise_counts.R $DiffAnalysis_WKD/TAPPAS_INPUT/ONT_unfiltered_expression.txt $RAWDIR/exp_factors.txt $DiffAnalysis_WKD/All/Unfiltered/SQANTI3/ONTTargeted_unfiltered_talon_classification.txt $DiffAnalysis_WKD/TAPPAS_INPUT/ONT_retained_list.txt $DiffAnalysis_WKD/TAPPAS_INPUT ONT

# File 3
grep -F -f $DiffAnalysis_WKD/TAPPAS_INPUT/ONT_retained_list.txt $DiffAnalysis_WKD/All/Unfiltered/SQANTI3/ONTTargeted_unfiltered_talon.gff3 > ONTTargeted_subset_talon.gff3

# File 4
cp $RAWDIR/ONT_phenotype.txt $DiffAnalysis_WKD/TAPPAS_INPUT

#exp=read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON//TAPPAS_INPUT/ONT_unfiltered_expression.txt", header = T)
#cols = c("K24", "L22", "M21" ,"O18", "O23", "O22", "P19", "T20", "Q21", "S18", "S23", "Q18", "Q17", "L18","T18")
#write.table(exp[,cols],"/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON//TAPPAS_INPUT/ONT_unfiltered_expression_missingQ23.txt", quote = F, sep = "\t")#

#retained_list = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON/TAPPAS_INPUT/ONT_retained_list.txt", header = T)
#expression = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON/TAPPAS_INPUT/ONT_normalised_expression.txt")

#final = subset(expression, rownames(expression) %in% retained_list$x)
#write.table(final, "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON/TAPPAS_INPUT/ONT_finalsubset_normalised_expression.txt", quote = F, sep = "\t")

################################################################################################
echo "#************************************* RNA-Seq Expression Matrix on Iso-Seq scaffold"
# all the samples RNASeq
SAMPLES_NAMES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)

# Prepare the Merged (Iso-Seq + ONT) fasta for aligment (with removed 3'ISM)
Rscript $GENERALFUNC/3ISM_remove_classification.R $mergedsqpath"_classification.txt" $mergedsqpath"_corrected.gtf" $mergedsqpath"_junctions.txt"  $DiffAnalysis_WKD/All/Merged/SQANTI3 IsoSeqONT
python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $mergedsqpath"_corrected.fasta" $DiffAnalysis_WKD/All/Merged/SQANTI3/IsoSeqONT_ISMrem.isoform.txt $mergedsqpath"_ISMrem.fasta"

# already merged from previous analysis (Diff_Whole.sh)
## 8) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
#mouse_merge_fastq $RNASeq_Filtered $DiffAnalysis_WKD/RNASeq_SQANTI3 AllRNASeq

## 9) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
source activate sqanti2
cd $DiffAnalysis_WKD/RNASeq_SQANTI3
kallisto index -i AllRNASeq_Kallisto.idx $mergedsqpath"_ISMrem.fasta" 2> AllRNASeq_Kallisto.index.log
SAMPLES_NAMES=(K24 L22 M20 O24 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 M18 O22 P20 Q18 S22 T20 K21 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 P23 S23 S17 T23)
for i in ${SAMPLES_NAMES[@]}; do run_kallisto_1sample $RNASeq_Filtered $i AllRNASeq_Kallisto.idx $DiffAnalysis_WKD/RNASeq_SQANTI3; done
##################################################################################################
#************************************* Run TappAS on Knight
scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON/TAPPAS_INPUT/* /mnt/data1/Szi/TAPPAS_ONT
scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON/TAPPAS_INPUT/ONT_finalsubset_normalised_expression.txt /mnt/data1/Szi/TAPPAS_ONT
/mnt/data1/Aaron/sw/jre1.8.0_181/bin/java -jar /mnt/data1/Szi/tappAS.jar

#### Tappas output
### tappAS projects
# Project.01394817487.tappas = Targeted Mouse + IsoSeq (All Trancripts for Input Expression)

scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.01394817487.tappas/* $DiffAnalysis_WKD/TAPPAS_OUTPUT/