#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=2_map_annotate_isoform.o1
#SBATCH --error=2_map_annotate_isoform.e1


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_isoseq_processing
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing

##-------------------------------------------------------------------------

# merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
merging_at_refine $WKD_ROOT/1_isoseq3/3_refine $WKD_ROOT/1_isoseq3/5_merged_cluster ${NAME} ${ALL_SAMPLE_NAMES[@]}
refine2fasta $WKD_ROOT/1_isoseq3/3_refine ${ALL_SAMPLE_NAMES[@]}

# align individual samples
# run_pbmm2align <output_name> <clustered_dir> <mapped_dir>
for i in ${ALL_SAMPLE_NAMES[@]}; do run_pbmm2align $i $WKD_ROOT/1_isoseq3/4_cluster $WKD_ROOT/2_post_isoseq3/6_minimap; done

# filter_alignment <name> <mapped_dir>
for i in ${ALL_SAMPLE_NAMES[@]}; do filter_alignment $i $WKD_ROOT/2_post_isoseq3/6_minimap; done

# run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
run_map_cupcakecollapse ${NAME} $WKD_ROOT/1_isoseq3/5_merged_cluster $WKD_ROOT/2_post_isoseq3/6_minimap $WKD_ROOT/2_post_isoseq3/7_tofu

# demux <name> <refine_dir> <cluster_report> <tofu_dir> 
demux ${NAME} $WKD_ROOT/1_isoseq3/3_refine $WKD_ROOT/1_isoseq3/5_merged_cluster/WholeIsoSeq.clustered.cluster_report.csv $WKD_ROOT/2_post_isoseq3/7_tofu


##-------------------------------------------------------------------------

# mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
mouse_merge_fastq ${RNASEQ_FILTERED_DIR} ${RNASEQ_MAPPED_DIR} ${NAME}

# run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto ${NAME} $WKD_ROOT/2_post_isoseq3/7_tofu/${NAME}.collapsed.fasta ${RNASEQ_MAPPED_DIR} $WKD_ROOT/2_post_isoseq3/8_kallisto


##-------------------------------------------------------------------------

# run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
run_sqanti3 ${NAME}.collapsed ${NAME}.collapsed.gff $WKD_ROOT/2_post_isoseq3/7_tofu ${RNASEQ_MAPPED_DIR} $WKD_ROOT/2_post_isoseq3/8_kallisto/${NAME}.mod.abundance.tsv $WKD_ROOT/2_post_isoseq3/7_tofu/${NAME}_fl_count.csv $WKD_ROOT/2_post_isoseq3/9_sqanti3 genome

# index for downstream alignment with rnaseq files
cd $WKD_ROOT/2_post_isoseq3/10_rnaseq2isoseq
export SQPATH=${WKD_ROOT}/2_post_isoseq3/9_sqanti3/${NAME}.collapsed
source activate sqanti2
kallisto index -i ${NAME}_Kallisto.idx ${SQPATH}.filtered.faa 2> ${NAME}_Kallisto.index.log

##-------------------------------------------------------------------------
# find human transgene sequences 

# create a file listing all the readID in cluster report (output = pre_cluster_read.csv)
mkdir ${WKD_ROOT}/2_post_isoseq3/11_transgene 
cd ${WKD_ROOT}/1_isoseq3/4_cluster/
for i in *cluster_report.csv*; do wc -l $i; done > pre_cluster_read.csv
sed -i 's/ \+/,/g' pre_cluster_read.csv 
sed -i '1 i num_reads,file' pre_cluster_read.csv
mv pre_cluster_read.csv ${WKD_ROOT}/2_post_isoseq3/11_transgene 


# grep transgene sequences in clustered .fasta
name=(hmapt1 hmapt2 mmapt1)
seq=($hMAPT_1 $hMAPT_2 $mMAPT_1)
for i in {0..2}; do
  echo ${name[$i]}
  echo ${seq[$i]}
  search_fasta_by_sequence.py --fasta=${WKD_ROOT}/1_isoseq3/4_cluster --i=hq.fa \
    --seq=${seq[$i]} -o=${name[$i]} \
    -d=${WKD_ROOT}/2_post_isoseq3/11_transgene 
done