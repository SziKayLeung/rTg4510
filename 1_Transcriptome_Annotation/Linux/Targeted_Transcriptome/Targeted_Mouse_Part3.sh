#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Targeted_Mouse_Part3.o
#SBATCH --error=Targeted_Mouse_Part3.e

# 03/06/2021: Same script as Targeted_Mouse_Part2.sh but only work with the same samples as Whole Transcriptome for fair comparisons

#************************************* DEFINE GLOBAL VARIABLES
# File directories
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/SUBSET
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/RNASeq
WHOLE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME

#cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER
#cd $PostIsoseq3_WKD; mkdir mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER Whole_vs_Targeted
#cd $RNASeq_WKD; mkdir MAPPED

# For Pooled Targeted
### Important order of BAM files is the same as the sample names
BATCH_NAMES=(Targeted_Seq_1 Targeted_Seq_2 Targeted_Seq_3a Targeted_Seq_3b)
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome
cat $RAWDIR/Isoseq_Targeted_MouseRaw.txt
BAM_FILES=(`cat $RAWDIR/Isoseq_Targeted_MouseRaw.txt | egrep -v "^\s*(#|$)"`)

# Other input files and directory
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019

# For Demultiplexing Samples
Pooled_Samples_Targeted1=(K19 K23 K21 K18 K20 K17)
Pooled_Samples_Targeted2=(S19 K24 L22 M21 O18 O23 O22 P19 T20)
Pooled_Samples_Targeted3a=(Q20A Q21A S18A S23A Q18A Q17A L18A Q23A T18A) # failed 3rd run
Pooled_Samples_Targeted3b=(Q20 Q21 S18 S23 Q18 Q17 L18 Q23 T18) # successful 3rd run
Barcoded_Targeted1_config_file=$RAWDIR/Barcode_Configs/Isoseq_Mouse_Targeted1_barcode.config
Barcoded_Targeted2_config_file=$RAWDIR/Barcode_Configs/Isoseq_Mouse_Targeted2_barcode.config
Barcoded_Targeted3_config_file=$RAWDIR/Barcode_Configs/Isoseq_Mouse_Targeted3_barcode.config
ALL_SAMPLES_NAMES=(K19 K23 K21 K18 K20 K17 S19 K24 L22 M21 O18 O23 O22 P19 T20 Q20 Q21 S18 S23 Q18 Q17 L18 Q23 T18)
ALL_TG4510_SAMPLES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)


# sourcing functions script
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome
source $FUNCTIONS/Targeted_Mouse_Functions.sh

module load Miniconda2/4.3.21
################################################################################################
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch (Mouse/Targeted_Mouse_Part1.sh) Functions 1 - 3, 6 for all batched runs"
echo "Already processed in individual samples for full set (Mouse/Targeted_Mouse_Part2.sh) Functions 4"

################################################################################################
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
## 5) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# Targeted_Seq_3b = complete run of Batch 3
merging_at_refine $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER SubsetAllMouseTargeted Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24

## 7) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
run_map_cupcakecollapse SubsetAllMouseTargeted $Isoseq3_WKD/MERGED_CLUSTER $PostIsoseq3_WKD/MAP $PostIsoseq3_WKD/TOFU

## 8) demux_targeted <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
demux_targeted $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER/SubsetAllMouseTargeted.clustered.cluster_report.csv $PostIsoseq3_WKD/TOFU/SubsetAllMouseTargeted.collapsed.read_stat.txt $PostIsoseq3_WKD/TOFU/SubsetAllMouseTargeted.Demultiplexed_Abundance.txt

################################################################################################
echo "#************************************* RNAseq [Function 9, 10]"
echo "Already processed in individual samples for rnaseq alignment for full set (Mouse/Targeted_Mouse_Part2.sh) Functions 9,10"

################################################################################################
echo "#************************************* RNASeq & IsoSeq [Function 11]"
## 11) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto AllMouseTargeted $PostIsoseq3_WKD/TOFU/SubsetAllMouseTargeted.collapsed.rep.fa $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO

################################################################################################
echo "#************************************* SQANTI2 [Function 12]"
## 12) run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
run_sqanti2 SubsetAllMouseTargeted.collapsed SubsetAllMouseTargeted.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/AllMouseTargeted.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/SubsetAllMouseTargeted.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2 genome

################################################################################################
echo "#************************************* TAMA filter [Function 13,14]"
## 13) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/SubsetAllMouseTargeted.collapsed_classification.filtered_lite.gtf SubsetAllMouseTargeted $PostIsoseq3_WKD/TAMA

## 14) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
sqname=SubsetAllMouseTargeted.collapsed_classification.filtered_lite
TAMA_sqanti_filter $PostIsoseq3_WKD/TAMA/SubsetAllMouseTargeted.bed $PostIsoseq3_WKD/SQANTI2 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" SubsetAllMouseTargeted $PostIsoseq3_WKD/SQANTI_TAMA_FILTER

################################################################################################
echo "#************************************* Whole vs Targeted Transcriptome"
## 18) convertgtf2bed <name> <working_directory> <input_sqanti_tama_directory>
convertgtf2bed WholeIsoSeq_sqantitamafiltered.classification $PostIsoseq3_WKD/Whole_vs_Targeted $WHOLE
convertgtf2bed SubsetAllMouseTargeted_sqantitamafiltered.classification $PostIsoseq3_WKD/Whole_vs_Targeted $PostIsoseq3_WKD/SQANTI_TAMA_FILTER

## 19) run_tamamerge <wkd> <wholebed> <targetedbed> <outputname>
run_tamamerge $PostIsoseq3_WKD/Whole_vs_Targeted $PostIsoseq3_WKD/Whole_vs_Targeted/WholeIsoSeq_sqantitamafiltered.classification.final.tama_mod.bed12 $PostIsoseq3_WKD/Whole_vs_Targeted/SubsetAllMouseTargeted_sqantitamafiltered.classification.final.tama_mod.bed12 merged_whole_targeted