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
#SBATCH --output=Differential_Analysis_All_Prepare.o
#SBATCH --error=Differential_Analysis_All_Prepare.e

# 25/01/2022: Strategy for differential analysis of whole transcriptome data

#### Strategy #####
# Long-read Iso-Seq data: SQANTI filtered dataset but with no RNA-Seq as input; remove only 3'ISM and counts (assumed partial products)
# Align short-read RNA-Seq data to long-read Iso-Seq data

#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis_Final
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq

mkdir -p $DiffAnalysis_WKD
cd $DiffAnalysis_WKD; mkdir -p SQANTI3 RNASeq_SQANTI3 TAPPAS_INPUT AS
cd $DiffAnalysis_WKD/TAPPAS_INPUT; mkdir -p IsoSeq_Expression RNASeq_Expression

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
R_FUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Rscripts
DIFF_FUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome

module load Miniconda2/4.3.21
source activate sqanti2_py3

dataset=WholeIsoSeq
sqname=$dataset".collapsed_classification.filtered_lite"
sqpath=$DiffAnalysis_WKD/SQANTI3/$sqname

################################################################################################
echo "#************************************* SQANTI3 [Function 1]"
## 1) run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
SAMPLES_NAMES=(Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)
run_sqanti3 $dataset.collapsed $dataset.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/$dataset".mod.abundance.tsv" $PostIsoseq3_WKD/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt $DiffAnalysis_WKD/SQANTI3 nornaseq

# Rscript script.R <input.classfile> <input.gtf> <output.dir> <prefix>
source activate sqanti2_py3
Rscript $GENERALFUNC/3ISM_remove_classification.R $sqpath"_classification.txt" $sqpath".gtf" $sqpath"_junctions.txt"  $DiffAnalysis_WKD/SQANTI3 $dataset
python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $sqpath".fasta" $DiffAnalysis_WKD/SQANTI3/$dataset"_ISMrem.isoform.txt" $sqpath"_ISMrem.fasta"

################################################################################################
echo "#************************************* RNA-Seq Expression Matrix on Iso-Seq scaffold"
# all the samples RNASeq
SAMPLES_NAMES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)

# already merged from previous analysis (Diff_Whole.sh)
## 8) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
#mouse_merge_fastq $RNASeq_Filtered $DiffAnalysis_WKD/RNASeq_SQANTI3 AllRNASeq

## 9) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# individual all the RNASeq samples
# first index the fasta file (note this made in run_kallisto but recreating for ease)
# later align all 59 RNASeq samples separately using Kallisto to IsoSeq scaffold (Diff_WholeRNASeq.sh)
source activate sqanti2
cd $DiffAnalysis_WKD/RNASeq_SQANTI3
kallisto index -i AllRNASeq_Kallisto.idx $sqpath"_ISMrem.fasta" 2> AllRNASeq_Kallisto.index.log
bash $DIFF_FUNC/Differential_Analysis_RNASeq_Prepare.sh
wait

##################################################################################################
echo "#************************************* Alternative Splicing Events for Groups"
source activate sqanti2_py3
Rscript $R_FUNC/Counts_Groups.R $sqpath"_ISM_rem.txt" $DiffAnalysis_WKD/AS

for group in WT_2mos TG_2mos WT_8mos TG_8mos WT TG; do
  echo "Processing: $group"
  # Rscript .R <retained_id> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_junc_txt> <output_prefix_name> <output_dir>
  Rscript $GENERALFUNC/sqanti_classgtfsubset.R $DiffAnalysis_WKD/AS/$group"_counts.txt" $DiffAnalysis_WKD/SQANTI3 $sqname"_ISMrem.classification" $sqname"_ISMrem.classification.gtf" $sqname"_ISMrem.junction.txt" $group $DiffAnalysis_WKD/AS

  # run_suppa2 <input_gtf> <input_class> <output_dir> <output_name>
  run_suppa2 $DiffAnalysis_WKD/AS/$group"_sqantisubset.classification.gtf" $DiffAnalysis_WKD/AS/$group"_sqantisubset.classification.txt" $DiffAnalysis_WKD/AS $group
done

