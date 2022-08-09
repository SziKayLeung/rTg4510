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
#SBATCH --output=Trial_RNASeqAlignment_Part2.o
#SBATCH --error=Trial_RNASeqAlignment_Part2.e

# 20/09/2021: Align 1 sample of RNA-Seq reads using different strategies

DiffAnalysis_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis
DIFF_FUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis
WKD=$DiffAnalysis_WKD/RNASeq_SQANTI3/TESTING
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/RNASeq
Whole_Targeted=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/DiffAnalysis/Whole_Targeted
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered

# sourcing functions script and input directories
module load Miniconda2/4.3.21
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/Whole_Transcriptome
source $FUNCTIONS/Diff_Whole_Functions.sh

source activate sqanti2_py3

prepare_alignment(){
  dataset=AllMouseTargeted
  sqname=$dataset".collapsed_classification.filtered_lite"

  # Script <sqanti_filered_inputfile> <tama_filtered_inputfile> <output_expression_file> <output_finalisoform_file> <type == "Targeted/Whole">
  Rscript $DIFF_FUNC/Sqanti_Collapsed_Counts.R $DiffAnalysis_WKD/SQANTI3/$sqname"_classification.txt" NA $WKD/AllReadsCollapsed $WKD/AllReadsCollapsed_Isoforms Whole

  Rscript $DIFF_FUNC/Sqanti_Collapsed_Counts.R $DiffAnalysis_WKD/SQANTI3/$sqname"_classification.txt" NA $WKD/TargetReadsCollapsed.txt $WKD/TargetReadsCollapsed_Isoforms Targeted

  for type in AllReadsCollapsed TargetReadsCollapsed; do
    echo $type
    length_name=$type"_FSMLength"
    Exp_name=$type"_FSMExp"

    #cd $WKD; mkdir $length_name $Exp_name

    collapse_filter $WKD/$type"_Isoforms_FSMbyexpression.txt" $DiffAnalysis_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" $dataset $WKD/$Exp_name

    collapse_filter $WKD/$type"_Isoforms_FSMbylength.txt" $DiffAnalysis_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" $dataset $WKD/$length_name
  done
}


kallisto_run(){
  index_file=$1
  input_fasta=$2
  input_gtf=$3
  kallisto index -i $1"_Kallisto.idx" $input_fasta 2> $1"_Kallisto.idx.log"
  kallisto quant -i $1"_Kallisto.idx" -rf-stranded $RNASeq_Filtered/Tg4510_filtered/K24/K24_S57_R1_001.fastq.filtered $RNASeq_Filtered/Tg4510_filtered/K24/K24_S57_R2_001.fastq.filtered -o K24_$1 2> $1"_Kallisto.quant.log"
  #kallisto quant -i $1"_Kallisto.idx" -b 30 -o K24_$1 --genomebam --gtf $input_gtf $RNASeq_Filtered/Tg4510_filtered/K24/K24_S57_R1_001.fastq.filtered $RNASeq_Filtered/Tg4510_filtered/K24/K24_S57_R2_001.fastq.filtered
}

kallisto_bootstrap_run(){
  index_file=$1
  input_fasta=$2
  input_gtf=$3
  kallisto index -i $1"_Kallisto.idx" $input_fasta 2> $1"_Kallisto.idx.log"
  kallisto quant -i $1"_Kallisto.idx" -b 200 -rf-stranded $RNASeq_Filtered/Tg4510_filtered/K24/K24_S57_R1_001.fastq.filtered $RNASeq_Filtered/Tg4510_filtered/K24/K24_S57_R2_001.fastq.filtered -o K24_$1 2> $1"_Kallisto.quant.log"
  #kallisto quant -i $1"_Kallisto.idx" -b 30 -o K24_$1 --genomebam --gtf $input_gtf $RNASeq_Filtered/Tg4510_filtered/K24/K24_S57_R1_001.fastq.filtered $RNASeq_Filtered/Tg4510_filtered/K24/K24_S57_R2_001.fastq.filtered
}

psuedobam_index(){
  wkd=$1  
  cd $wkd
  samtools index pseudoalignments.bam -csi -m 14 
}

# prepare_alignment

# 1. All Targeted Reads; Not Collapsed
# 2. All Targeted Reads; FSM Only (Chosen by longest transcript)
# 3. On target genes only; FSM Only (Chosen by longest transcript)
# 4. All Targeted Reads; FSM Only (Chosen by most abundant transcript)
# 5. On target genes only; FSM Only (Chosen by most abundant transcript)
# 6. Merged Whole & Targeted; FSM Only (Chosen by longest transcript)

fasta=AllMouseTargeted_sqantisubset_classification.fasta
gtf=AllMouseTargeted_sqantisubset.final.classification.gtf

cd $WKD
kallisto_run AllReadsNotCollapsed $DiffAnalysis_WKD/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite.fasta $DiffAnalysis_WKD/SQANTI3/AllMouseTargeted.collapsed_classification.filtered_lite.gtf

kallisto_run AllReadsCollapsed_FSMLength $WKD/AllReadsCollapsed_FSMLength/$fasta $WKD/AllReadsCollapsed_FSMLength/$gtf
kallisto_run TargetedReadsCollapsed_FSMLength $WKD/TargetReadsCollapsed_FSMLength/$fasta $WKD/TargetReadsCollapsed_FSMLength/$gtf
kallisto_run AllReadsCollapsed_FSMExp $WKD/AllReadsCollapsed_FSMExp/$fasta $WKD/AllReadsCollapsed_FSMExp/$gtf
kallisto_run TargetedReadsCollapsed_FSMExp $WKD/AllReadsCollapsed_FSMExp/$fasta $WKD/AllReadsCollapsed_FSMExp/$gtf
kallisto_run Merged_WholeTargeted $Whole_Targeted/merged_whole_targeted_corrected.fasta $Whole_Targeted/merged_whole_targeted_corrected.gtf

pseudoalignments bam 
dataset=(AllReadsNotCollapsed AllReadsCollapsed_FSMLength TargetedReadsCollapsed_FSMLength AllReadsCollapsed_FSMExp TargetedReadsCollapsed_FSMExp Merged_WholeTargeted) 
for d in ${dataset[@]}; do 
  echo "Processing psuedalignments.bam in $WKD/K24_$d"
  cd $WKD/"K24_"$d 
  samtools index pseudoalignments.bam -c -m 14  
done

kallisto_bootstrap_run AllReadsCollapsed_FSMExp_BTS $WKD/AllReadsCollapsed_FSMExp/$fasta $WKD/AllReadsCollapsed_FSMExp/$gtf

##### ERCC 
#ERCC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/ERCC

#kallisto_run ERCCReadsNotCollapsed $ERCC/WholeIsoSeq.collapsed_classification.filtered_lite.fasta $ERCC/WholeIsoSeq.collapsed_classification.filtered_lite.gtf
#kallisto_run ERCCReadsCollapsed $ERCC/WholeIsoSeq_sqantifiltered_tamafiltered_classification.fasta $ERCC/WholeIsoSeq_sqantitamafiltered.classification.gtf
