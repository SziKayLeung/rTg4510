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
#SBATCH --array=0-11 # 12 samples
#SBATCH --output=Individual_Part1-%A_%a.o
#SBATCH --error=Individual_Part1-%A_%a.e

# 13/12/2020: Running all 12 Tg4510 samples separately in batches from CCS to cupcake
# 14/01/2021: Rerun mapping but with stats output
# 15/01/2021: Mapping of individual samples to ERCC for mapping quality

#************************************* DEFINE GLOBAL VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome
#cd $Isoseq3_WKD
#mkdir CCS LIMA REFINE CLUSTER MAPPING TOFU ERCC_MAPPING
CCS=$Isoseq3_WKD/CCS
LIMA=$Isoseq3_WKD/LIMA
REFINE=$Isoseq3_WKD/REFINE
CLUSTER=$Isoseq3_WKD/CLUSTER
MAPPING=$Isoseq3_WKD/MAPPING
TOFU=$Isoseq3_WKD/TOFU
ERCC_MAPPING=$Isoseq3_WKD/ERCC_MAPPING

# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
SAMPLES_NAMES=(Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)
SAMPLE=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
cd $RAWDIR
cat Isoseq_Whole_MouseRaw.txt
# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
BAM_FILES=(`cat "Isoseq_Whole_MouseRaw.txt" | egrep -v "^\s*(#|$)"`)
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}

#************************************* IsoSeq - ALL samples separately
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
# Isoseq3.4.0
    # run_CCS_batch <input_ccs_bam> <prefix_output_name> <Output_directory>
    # run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
    # run_REFINE $Sample $Input_LIMA_directory $Output_directory
    # run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
#run_CCS_batch ${BAM_FILE} ${SAMPLE} $CCS
#run_LIMA ${SAMPLE} $CCS $LIMA "no_multiplex"
#run_REFINE ${SAMPLE} $LIMA $REFINE
#run_CLUSTER ${SAMPLE} $REFINE $CLUSTER

#python $FUNCTIONS/IsoSeq_QC/CCS.py $CCS "" All
#python $FUNCTIONS/IsoSeq_QC/LIMA.py $LIMA "" All

#************************************* Post-IsoSeq - ALL samples separately
source $FUNCTIONS/Post_Isoseq3_Function.sh
# convert_fa2fq <file_name> <input_dir>
# run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
# tofu <prefix_sample> <cluster_dir> <input_dir> <output_dir>
#convert_fa2fq ${SAMPLE}.clustered.hq.fasta $CLUSTER
run_minimap2 ${SAMPLE} $CLUSTER mm10 $MAPPING
run_minimap2 ${SAMPLE} $CLUSTER ERCC $ERCC_MAPPING

tofu_re(){
    cd $4; mkdir $1
    cd $1
    source activate cupcake
    echo "Processing Sample $1 for TOFU"
    # Collapse
    collapse_isoforms_by_sam.py --input $2/$1.clustered.hq.fastq --fq -s $3/$1.clustered.hq.fastq.sorted.sam --dun-merge-5-shorter -o Sample &> $1.collapse.log

    # Create Abundance Script of full-length transcripts
    get_abundance_post_collapse.py Sample.collapsed $2/$1.clustered.cluster_report.csv 2> $1.abundance.log

    # Remove degraded isoforms (default setting)
    filter_away_subset.py Sample.collapsed 2> $1.filter.log

    source deactivate

    source activate sqanti2_py3
    # convert rep.fq to rep.fa for SQANTI2 input
    seqtk seq -a Sample.collapsed.filtered.rep.fq > Sample.collapsed.filtered.rep.fa
    echo "Processing Sample $1 for TOFU successful"
}

#tofu_re ${SAMPLE} $CLUSTER $MAPPING $TOFU
