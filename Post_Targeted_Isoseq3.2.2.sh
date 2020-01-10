#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrcq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 10/01/2019: Created script to run Post_Isoseq3.2.2 on Mouse Targeted Run 1 (TargetedSeq1)
    # Samples K20, K18, K23, K21, K19, K17 

sl693=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693
FUNCTIONS=$sl693/Scripts/general


#************************************* DEFINE VARIABLES
TS1_dir=$sl693/Targeted/Targeted_Seq_1
SAMPLES_NAMES=(TargetedSeq1)
samples=(K18 K23 K20 K21 K19 K17)
POLISHED=$TS1_dir/CLUSTER
MAPPING=$TS1_dir/MAPPING
TOFU=$TS1_dir/TOFU
SQANTI2_output_dir=$TS1_dir/SQANTI2
SQANTI2_filter_output_dir=$TS1_dir/SQANTI2
RNASEQ=$sl693/RNASeq/all_filtered
STAR=$TS1_dir/RNASEQ

#************************************* Run script
source $FUNCTIONS/RNASeq/STAR_Functions.sh

# Define function for STAR merge using parameters for junction reads
# run_star_merge <J20/Tg4510_input_directory> <WKD> <sample_name_output>
# DEFINE samples=(sample_1, sample_2, sample_n) in working script
run_star_merge $RNASEQ $TS1_dir/RNASEQ TargetedSeq1

source $FUNCTIONS/Post_IsoSeq/Post_Isoseq3_Functions.sh

# convert_fa2fq <sample_name> <input_dir>
# run_minimap2 $sample_1 $sample_2
# tofu $sample_1 $sample_2

for i in ${SAMPLES_NAMES[@]}; do
    convert_fa2fq $i $POLISHED
    run_minimap2 $i 
    tofu $i
    run_sqanti2_QC $i
    run_sqanti2_Filter $i
done

PROBES=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/Expected_Expression/mouse
python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/targeted/calc_probe_hit_from_sam.py \
         $PROBES/FINAL_MOUSE.bed $POLISHED/TargetedSeq1.polished.hq.fasta $MAPPING/TargetedSeq1.polished.hq.fastq.sam \
         --start_base 0 --end_base 0 \
         -o TargetedSeq1.fasta.sam.probe_hit.txt
