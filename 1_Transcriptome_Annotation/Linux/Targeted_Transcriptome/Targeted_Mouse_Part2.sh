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
#SBATCH --output=Targeted_Mouse_Part2.o
#SBATCH --error=Targeted_Mouse_Part2.e

#************************************* DEFINE GLOBAL VARIABLES
# File directories
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/RNASeq

cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER
cd $PostIsoseq3_WKD; mkdir mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER
cd $RNASeq_WKD; mkdir MAPPED

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

# sourcing functions script
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome
source $FUNCTIONS/Targeted_Mouse_Functions.sh

module load Miniconda2/4.3.21
################################################################################################
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch (Mouse/Targeted_Mouse_Part1.sh) Functions 1 - 3, 6 for all batched runs"

## 4) run_targeted_REFINE ${BATCH} $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
for i in ${Pooled_Samples_Targeted1[@]}; do run_targeted_REFINE $i $Barcoded_Targeted1_config_file Targeted_Seq_1 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done
for i in ${Pooled_Samples_Targeted2[@]}; do run_targeted_REFINE $i $Barcoded_Targeted2_config_file Targeted_Seq_2 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done
for i in ${Pooled_Samples_Targeted3a[@]}; do run_targeted_REFINE $i $Barcoded_Targeted3_config_file Targeted_Seq_3a $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done
for i in ${Pooled_Samples_Targeted3b[@]}; do run_targeted_REFINE $i $Barcoded_Targeted3_config_file Targeted_Seq_3b $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done

# 6) run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
for i in ${ALL_SAMPLES_NAMES[@]}; do run_CLUSTER $i $Isoseq3_WKD/REFINE $Isoseq3_WKD/CLUSTER; done

################################################################################################
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
## 5) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# Targeted_Seq_3b = complete run of Batch 3
merging_at_refine $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER AllMouseTargeted Targeted_Seq_1 Targeted_Seq_2 Targeted_Seq_3b

## 7) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
run_map_cupcakecollapse AllMouseTargeted $Isoseq3_WKD/MERGED_CLUSTER $PostIsoseq3_WKD/MAP $PostIsoseq3_WKD/TOFU

## 8) demux_targeted <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
demux $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER/AllMouseTargeted.clustered.cluster_report.csv $PostIsoseq3_WKD/TOFU/AllMouseTargeted.collapsed.read_stat.txt $PostIsoseq3_WKD/TOFU/AllMouseTargeted.Demultiplexed_Abundance.txt

################################################################################################
echo "#************************************* RNAseq [Function 9, 10]"
## 9) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
for i in ${ALL_SAMPLES_NAMES[@]}; do run_star $i $RNASeq_Filtered $RNASeq_WKD/MAPPED; done
cd $RNASeq_WKD/MAPPED; for i in ${ALL_SAMPLES_NAMES[@]}; do echo "################## $i"; cat $i/$i.Log.final.out;done > Mapped_Results.txt

## 10) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
mouse_merge_fastq $RNASeq_Filtered $RNASeq_WKD/MAPPED AllMouseTargeted

################################################################################################
echo "#************************************* RNASeq & IsoSeq [Function 11]"
## 11) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto AllMouseTargeted $PostIsoseq3_WKD/TOFU/AllMouseTargeted.collapsed.rep.fa $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO

################################################################################################
echo "#************************************* SQANTI2 [Function 12]"
## 12) run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
run_sqanti2 AllMouseTargeted.collapsed AllMouseTargeted.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/AllMouseTargeted.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/AllMouseTargeted.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2/GENOME genome

################################################################################################
echo "#************************************* TAMA filter [Function 13,14]"
## 13) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/GENOME/AllMouseTargeted.collapsed_classification.filtered_lite.gtfAllMouseTargeted $PostIsoseq3_WKD/TAMA/GENOME

## 14) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
sqname=AllMouseTargeted.collapsed_classification.filtered_lite
TAMA_sqanti_filter $PostIsoseq3_WKD/TAMA/GENOME/AllMouseTargeted.bed $PostIsoseq3_WKD/SQANTI2/GENOME $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt"AllMouseTargeted $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/GENOME

################################################################################################
echo "#************************************* QC [Function 15]"
# parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
parse_stats_per_sample $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA AllMouseTargeted
