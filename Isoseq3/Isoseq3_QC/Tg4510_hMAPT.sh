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

# 22/09/2020: Align all individual Tg4510 samples to hg38 genome for MAPT validation

#************************************* DEFINE GLOBAL VARIABLES
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general
POLISHED=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER
MAPT_VALIDATION=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER/hMAPT
TO_HUMAN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/TO_HUMAN
MAPPING=$TO_HUMAN/MAPPING
TOFU=$TO_HUMAN/TOFU
SQANTI2_output_dir=$TO_HUMAN/SQANTI2_v7
STAR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/MAPPED/Individual
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered

#************************************* METHOD1: Occurence of human MAPT transgene in CCS reads 
hMAPT_1=TGGTTAATCACTTAACCTGCTTTTGTCACTCGGCTTTGGCTCGGGACTTCAAAATCAGTGATGGGAGTAAGAGCAAATTTCATCTTTCCAAATTGATGGGTGGGCTAGTAATAAAATATTTAAAAAAAAACATTCAAAAACATGGCCACATCCAACATTTCCTCAGGCAATTCCTTTTGATTCTTTTTTCTTCCCCCTCCATGTA
hMAPT_2=AAAATCAGTGATGGGAGTAAGAGCAAATTTCATCTTTCCAAATTGATGGGTGGGCTAGTAATAAAATATTTAAAAAAAAACATTCAAAAACATGGCCACATCCAACATTTCCTCAGGCAATTCCTTTTGATTCTTTTTTCTTCCCCCTCCATGTAGAAGAGGGAGAAGGAGAGGCTCTGAAAGCTGCTTCTGGGGGATTT

cd $POLISHED
for i in *fastq; do 
	echo "Processing with $i"
	sample=$(basename "$i" | cut -d "." -f 1)
	
	grep -B1 $hMAPT_1 $i > $MAPT_VALIDATION/$sample.hMAPT1.polished.hq.fastq
	grep "^@" $MAPT_VALIDATION/$sample.hMAPT1.polished.hq.fastq > $MAPT_VALIDATION/$sample.hMAPT1.header
	
	grep -B1 $hMAPT_2 $i > $MAPT_VALIDATION/$sample.hMAPT2.polished.hq.fastq
	grep "^@" $MAPT_VALIDATION/$sample.hMAPT2.polished.hq.fastq > $MAPT_VALIDATION/$sample.hMAPT2.header
done 

#************************************* METHOD12:  Gene expression of MAPT in final mouse transcriptome 
SAMPLES_NAMES=(Q21 O18 C21 E18 C20 B21 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)
source $FUNCTIONS/Post_IsoSeq/Post_Isoseq3_Functions.sh
for i in ${SAMPLES_NAMES[@]}; do 
    convert_fa2fq $i $POLISHED
    run_minimap2 $i $POLISHED hg38
    tofu $i $POLISHED
	
	# run_sqanti2_QC <prefix_sample> <input_tofu_dir> <coverage/genome=mm10_rnqaseq/mm10/hg38_gencode/hg38_chess> <input_kallisto_file> <output_dir> <input_rnaseq_dir> 
    run_sqanti2_QC $i $TOFU hg38_gencode NA $SQANTI2_output_dir NA 
    run_sqanti2_Filter $i 
done
