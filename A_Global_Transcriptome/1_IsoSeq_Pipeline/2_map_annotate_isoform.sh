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
#SBATCH --output=Targeted_ADBDR_Part2c.o
#SBATCH --error=Targeted_ADBDR_Part2c.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/A_Global_Transcriptome
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing

##-------------------------------------------------------------------------

# merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
merging_at_refine $WKD_ROOT/1_isoseq3/3_refine $WKD_ROOT/1_isoseq3/5_merged_cluster ${NAME} ${ALL_SAMPLE_NAMES[@]}
refine2fasta $WKD_ROOT/1_isoseq3/3_refine ${ALL_SAMPLE_NAMES[@]}

# run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
run_map_cupcakecollapse ${NAME} $WKD_ROOT/1_isoseq3/5_merged_cluster $WKD_ROOT/2_post_isoseq3/6_minimap $WKD_ROOT/2_post_isoseq3/7_tofu

# demux <input path read.stat file> <input path of samples file> <path of output>
demux $WKD_ROOT/2_post_isoseq3/7_tofu/${NAME}.collapsed.read_stat.txt ${SAMPLE_CONFIG} $WKD_ROOT/2_post_isoseq3/7_tofu/${NAME}.Demultiplexed_Abundance.txt


##-------------------------------------------------------------------------

# mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
mouse_merge_fastq ${RNASEQ_FILTERED_DIR} ${RNASEQ_MAPPED_DIR} ${NAME}

# run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto ${NAME} $WKD_ROOT/2_post_isoseq3/7_tofu/${NAME}.collapsed.rep.fa ${RNASEQ_MAPPED_DIR} $WKD_ROOT/2_post_isoseq3/8_kallisto


##-------------------------------------------------------------------------

# run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
run_sqanti3 ${NAME}.collapsed ${NAME}.collapsed.gff $WKD_ROOT/2_post_isoseq3/7_tofu ${RNASEQ_MAPPED_DIR} $WKD_ROOT/2_post_isoseq3/8_kallisto/${NAME}.mod.abundance.tsv $WKD_ROOT/2_post_isoseq3/7_tofu/${NAME}.Demultiplexed_Abundance.txt $DiffAnalysis_WKD/SQANTI3 genome

# Rscript script.R <input.classfile> <input.gtf> <output.dir> <prefix>
source activate sqanti2_py3
Rscript ${ISMREMOVE} ${SQPATH}_classification.txt ${SQPATH}.gtf ${SQPATH}_junctions.txt $WKD_ROOT/2_post_isoseq3/9_sqanti3 ${NAME}
python ${TAMASUBSETFASTA} ${SQPATH}.fasta $WKD_ROOT/2_post_isoseq3/9_sqanti3/${NAME}_ISMrem.isoform.txt ${SQPATH}_ISMrem.fasta

# index for downstream alignment with rnaseq files
cd $WKD_ROOT/2_post_isoseq3/10_rnaseq2isoseq
kallisto index -i AllRNASeq_Kallisto.idx ${SQPATH}_ISMrem.fasta 2> AllRNASeq_Kallisto.index.log


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