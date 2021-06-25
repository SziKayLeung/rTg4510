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
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/RNASeq
WholevsTargeted=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Whole_vs_Targeted
WHOLE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME

#cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER
#cd $PostIsoseq3_WKD; mkdir mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER SQANTI2_allrnaseq HUMANMAPT TAPPAS ALLRNASEQ
#cd $RNASeq_WKD; mkdir MAPPED

# For Pooled Targeted
### Important order of BAM files is the same as the sample names
BATCH_NAMES=(Targeted_Seq_1 Targeted_Seq_2 Targeted_Seq_3a Targeted_Seq_3b)
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome
Targeted_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/
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
SAMPLES_NAMES=(K19 K23 K21 K18 K20 K17 S19 K24 L22 M21 O18 O23 O22 P19 T20 Q20 Q21 S18 S23 Q18 Q17 L18 Q23 T18)
ALL_TG4510_SAMPLES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)


# sourcing functions script
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome
source $FUNCTIONS/Targeted_Mouse_Functions.sh

module load Miniconda2/4.3.21
################################################################################################
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch (Mouse/Targeted_Mouse_Part1.sh) Functions 1 - 3, 6 for all batched runs"

for i in ${BATCH_NAMES[@]};do run_LIMA $i $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA multiplex;done

## 4) run_targeted_REFINE ${BATCH} $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
for i in ${Pooled_Samples_Targeted1[@]}; do run_targeted_REFINE $i $Barcoded_Targeted1_config_file Targeted_Seq_1 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done
for i in ${Pooled_Samples_Targeted2[@]}; do run_targeted_REFINE $i $Barcoded_Targeted2_config_file Targeted_Seq_2 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done
for i in ${Pooled_Samples_Targeted3b[@]}; do run_targeted_REFINE $i $Barcoded_Targeted3_config_file Targeted_Seq_3b $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done

# 6) run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
for i in ${SAMPLES_NAMES[@]}; do run_CLUSTER $i $Isoseq3_WKD/REFINE $Isoseq3_WKD/CLUSTER; done

################################################################################################
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8, 9]"
## 5) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# Targeted_Seq_3b = complete run of Batch 3
merging_at_refine $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER AllMouseTargeted ${SAMPLES_NAMES[@]}

## 7) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
run_map_cupcakecollapse AllMouseTargeted $Isoseq3_WKD/MERGED_CLUSTER $PostIsoseq3_WKD/MAP $PostIsoseq3_WKD/TOFU

## 8) demux_targeted <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
demux_targeted $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER/AllMouseTargeted.clustered.cluster_report.csv $PostIsoseq3_WKD/TOFU/AllMouseTargeted.collapsed.read_stat.txt $PostIsoseq3_WKD/TOFU/AllMouseTargeted.Demultiplexed_Abundance.txt

## 9) on_target_rate <sample> <input_cluster_dir> <output_dir> <input_probe_bed.file>
# command line 
cd $PostIsoseq3_WKD/MAP; mkdir TARGETRATE
for i in ${SAMPLES_NAMES[@]};do on_target_rate $i $Isoseq3_WKD/CLUSTER $PostIsoseq3_WKD/MAP/TARGETRATE $RAWDIR/Probes/FINAL_MOUSE.bed; done

################################################################################################
echo "#************************************* RNAseq [Function 10,11]"
## 10) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
for i in ${SAMPLES_NAMES[@]}; do run_star $i $RNASeq_Filtered $RNASeq_WKD/MAPPED; done
cd $RNASeq_WKD/MAPPED; for i in ${SAMPLES_NAMES[@]}; do echo "################## $i"; cat $i/$i.Log.final.out;done > Mapped_Results.txt

# output: Targeted_Mouse_Part2Rnaseq.o; error: Targeted_Mouse_Part2Rnaseq.e
for i in ${ALL_TG4510_SAMPLES[@]}; do run_star $i $RNASeq_Filtered $RNASeq_WKD/MAPPED; done

## 11) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
SAMPLES_NAMES=(K19 K23 K21 K18 K20 K17 S19 K24 L22 M21 O18 O23 O22 P19 T20 Q20 Q21 S18 S23 Q18 Q17 L18 Q23 T18)
mouse_merge_fastq $RNASeq_Filtered $RNASeq_WKD/MAPPED AllMouseTargeted

################################################################################################
echo "#************************************* RNASeq & IsoSeq [Function 12]"
## 12) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto AllMouseTargeted $PostIsoseq3_WKD/TOFU/AllMouseTargeted.collapsed.rep.fa $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO

################################################################################################
echo "#************************************* SQANTI2 [Function 13]"
## 13) run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
run_sqanti2 AllMouseTargeted.collapsed AllMouseTargeted.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/AllMouseTargeted.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/AllMouseTargeted.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2 genome

################################################################################################
echo "#************************************* TAMA filter [Function 14,15]"
## 14) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/AllMouseTargeted.collapsed_classification.filtered_lite.gtf AllMouseTargeted $PostIsoseq3_WKD/TAMA

## 15) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
sqname=AllMouseTargeted.collapsed_classification.filtered_lite
TAMA_sqanti_filter $PostIsoseq3_WKD/TAMA/AllMouseTargeted.bed $PostIsoseq3_WKD/SQANTI2 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" AllMouseTargeted $PostIsoseq3_WKD/SQANTI_TAMA_FILTER

# TAMA remove fragments from prefiltered SQANTI datasets
cd $PostIsoseq3_WKD/SQANTI2; mkdir TAMAFILTERSQANTI
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/AllMouseTargeted.collapsed_corrected.gtf AllMouseTargeted $PostIsoseq3_WKD/SQANTI2/TAMAFILTERSQANTI
sqname=AllMouseTargeted.collapsed
TAMA_sqanti_filter $PostIsoseq3_WKD/SQANTI2/TAMAFILTERSQANTI/AllMouseTargeted.bed $PostIsoseq3_WKD/SQANTI2 $sqname"_classification.txt" $sqname"_corrected.gtf" $sqname"_corrected.fasta" $sqname"_junctions.txt" AllMouseTargeted $PostIsoseq3_WKD/SQANTI2/TAMAFILTERSQANTI

################################################################################################
echo "#************************************* QC [Function 16]"
# 16) parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
parse_stats_per_sample $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA AllMouseTargeted

################################################################################################
echo "#************************************* Repeat SQANTI2 with more RNASeq coverage [Function 11, 12, 13, 14, 15]"
SAMPLES_NAMES=(K24 L22 M20 O24 P22 Q20 S24 T22 K17 L21 M19 K23 P21 Q19 M21 T21 K22 L20 M18 O22 P20 Q18 S22 T20 K21 L19 M17 O21 P19 Q17 S21 T19 K20 L18 M24 O20 P18 Q24 S20 T18 K19 L17 M23 O19 P17 Q23 S19 T17 K18 L24 M22 O18 P24 Q22 S18 T24 O23 L23 Q21 O17 P23 S23 S17 T23)
mouse_merge_fastq $RNASeq_Filtered $RNASeq_WKD/MAPPED AllMouse
run_sqanti2 AllMouseTargeted.collapsed AllMouseTargeted.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/AllMouseTargeted.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/AllMouseTargeted.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2_allrnaseq genome
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2_allrnaseq/AllMouseTargeted.collapsed_classification.filtered_lite.gtf AllMouseTargeted $PostIsoseq3_WKD/SQANTI2_allrnaseq
sqname=AllMouseTargeted.collapsed_classification.filtered_lite
TAMA_sqanti_filter  $PostIsoseq3_WKD/SQANTI2_allrnaseq/AllMouseTargeted.bed $PostIsoseq3_WKD/SQANTI2_allrnaseq $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" AllMouseTargeted $PostIsoseq3_WKD/SQANTI2_allrnaseq

################################################################################################
echo "#************************************* Find human MAPT [Function 17]"
# 17) find_humanMAPT <cluster_dir> <output_dir> <merged_cluster.fa>
find_humanMAPT $Isoseq3_WKD/CLUSTER $PostIsoseq3_WKD/HUMANMAPT $Isoseq3_WKD/MERGED_CLUSTER/AllMouseTargeted.clustered.hq.fasta

################################################################################################
echo "#************************************* Whole vs Targeted Transcriptome [Function 18,19]"
# In Targeted_Part3.sh

################################################################################################
echo "#************************************* Differential Analysis - tappAS [Function 20,21]"
# 20) isoannolite_generate <input_dir> <input_gtf> <input_class> <input_junc> <species> <output_dir> <output_name>
# generate IsoAnnotLite output required for TAPPAS after sqanti (final files )
isoannolite_generate $PostIsoseq3_WKD/SQANTI_TAMA_FILTER AllMouseTargeted_sqantitamafiltered.classification.final.gtf AllMouseTargeted_sqantitamafiltered.classification.txt AllMouseTargeted_sqantitamafiltered.junction.txt Mouse $PostIsoseq3_WKD/TAPPAS AllMouseTargeted_tappasannot_from_SQANTI2.gff3

# 21) counts_subset_4tappas <input_class> <output_class> <type_genes>
counts_subset_4tappas $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/AllMouseTargeted_sqantitamafiltered.classification.txt  $PostIsoseq3_WKD/TAPPAS/AllMouseTargeted_sqantitamafiltered.expression.txt AD

#### tappas on knight
#cd /mnt/data1/Szi/TAPPAS_MouseTargeted
#scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/TAPPAS/* .
#scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetedMouse_PhenotypeTAPPAS.txt .
#scp -r sl693@login.isca.ex.ac.uk:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Post_IsoSeq/ALLRNASEQ/TargetedMouseRNASeq_sqantitamafiltered.expression.txt .
#/mnt/data1/Aaron/sw/jre1.8.0_181/bin/java -jar /mnt/data1/Szi/tappAS.jar

#### after running tappas on knight
#cd $PostIsoseq3_WKD/TAPPAS; mkdir Results TAPPAS_output
#cd $PostIsoseq3_WKD/TAPPAS/Results; mkdir RNASeq IsoSeq
#cd $PostIsoseq3_WKD/TAPPAS/TAPPAS_output; mkdir RNASeq IsoSeq
#scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/TAPPAS_MouseTargeted/*tsv* $PostIsoseq3_WKD/TAPPAS/Results/IsoSeq
#### note remove the comment from the header for each tsv
#scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.546554991.tappas/InputData/input_normalized_matrix.tsv $PostIsoseq3_WKD/TAPPAS/Results/IsoSeq
#scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.0160514833.tappas/InputData/input_normalized_matrix.tsv $PostIsoseq3_WKD/TAPPAS/Results/RNASeq
#scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/TAPPAS_MouseTargeted/TargetedMouse+RNASeq/*tsv* $PostIsoseq3_WKD/TAPPAS/Results/RNASeq

##### tappAS projects
## Project.546554991.tappas = Targeted Mouse + IsoSeq
## Project.0160514833.tappas = Targeted Mouse + RNASeq
#scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.546554991.tappas/Data/time_factors.txt $PostIsoseq3_WKD/TAPPAS/TAPPAS_output
#scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.546554991.tappas/Data/gene_matrix.tsv $PostIsoseq3_WKD/TAPPAS/TAPPAS_output
#scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.546554991.tappas/Data/transcript_matrix.tsv $PostIsoseq3_WKD/TAPPAS/TAPPAS_output

#for file in time_factors.txt gene_matrix.tsv transcript_matrix.tsv; do
#  echo "Copying $file from Project.546554991.tappas to TAPPAS_output/IsoSeq"
#  scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.546554991.tappas/Data/$file $PostIsoseq3_WKD/TAPPAS/TAPPAS_output/IsoSeq
#  echo "Copying $file from Project.0160514833.tappas to TAPPAS_output/RNASeq"
#  scp -r sLeung@knight.ex.ac.uk:/home/sLeung/tappasWorkspace/Projects/Project.0160514833.tappas/Data/$file $PostIsoseq3_WKD/TAPPAS/TAPPAS_output/RNASeq
#done
