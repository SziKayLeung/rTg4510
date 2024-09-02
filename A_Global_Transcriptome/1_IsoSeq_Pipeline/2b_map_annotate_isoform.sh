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
#SBATCH --output=2b_map_annotate_isoform.o2
#SBATCH --error=2b_map_annotate_isoform.e2

# J20 C20, C21, B21 and E18 Iso-Seq pipeline


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/A_Global_Transcriptome
LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
source $SC_ROOT/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/1_IsoSeq_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_isoseq_processing
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing

J20_ALL_SAMPLE_NAMES=(B21 C20 C21 E18)

##-------------------------------------------------------------------------

# merging_at_refine <input_flnc_bam_dir> <output_directory> <output_J20NAME> <samples.....>
#merging_at_refine $WKD_ROOT/1_isoseq3/3_refine $WKD_ROOT/1_isoseq3/5_merged_cluster ${J20NAME} ${J20_ALL_SAMPLE_NAMES[@]}
#refine2fasta $WKD_ROOT/1_isoseq3/3_refine ${J20_ALL_SAMPLE_J20NAMES[@]}

# align individual samples
# run_pbmm2align <output_J20NAME> <clustered_dir> <mapped_dir>
#for i in ${J20_ALL_SAMPLE_NAMES[@]}; do run_pbmm2align $i $WKD_ROOT/1_isoseq3/4_cluster $WKD_ROOT/2_post_isoseq3/6_minimap; done

# filter_alignment <J20NAME> <mapped_dir>
#for i in ${J20_ALL_SAMPLE_NAMES[@]}; do filter_alignment $i $WKD_ROOT/2_post_isoseq3/6_minimap; done

# run_map_cupcakecollapse <sample_prefix_input/output_J20NAME> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
#run_map_cupcakecollapse ${J20NAME} $WKD_ROOT/1_isoseq3/5_merged_cluster $WKD_ROOT/2_post_isoseq3/6_minimap $WKD_ROOT/2_post_isoseq3/7_tofu

# demux <J20NAME> <refine_dir> <cluster_report> <tofu_dir> 
#demux ${J20NAME} $WKD_ROOT/1_isoseq3/3_refine $WKD_ROOT/1_isoseq3/5_merged_cluster/${J20NAME}.clustered.cluster_report.csv $WKD_ROOT/2_post_isoseq3/7_tofu


##-------------------------------------------------------------------------

source activate sqanti2_py3
cd $WKD_ROOT/2_post_isoseq3/9_sqanti3
python ${SQANTI3_DIR}/sqanti3_qc.py -t 30 $WKD_ROOT/2_post_isoseq3/7_tofu/${J20NAME}.collapsed.gff ${GENOME_GTF} ${GENOME_FASTA} --CAGE_peak ${CAGE_PEAK} --polyA_motif_list ${POLYA} --genename --isoAnnotLite --report skip &> ${J20NAME}.collapsed.sqanti.qc.log
    