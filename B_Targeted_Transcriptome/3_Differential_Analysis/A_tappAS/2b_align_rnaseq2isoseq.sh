#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

## ---------------------------
## Purpose: Align RNA-Seq (core 96 samples) to Iso-Seq annotation for differential expression analysis
## ---------------------------

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/B_Targeted_Transcriptome
source $SC_ROOT/3_Differential_Analysis/rTg4510_differential.config
source $SC_ROOT/3_Differential_Analysis/01_source_functions.sh


##-------------------------------------------------------------------------
echo "#************************************* RNA-Seq Expression Matrix on Iso-Seq scaffold"

mkdir -p $HYBRID_DIFF_DIR/alignment

# index Iso-Seq fasta file
cd $HYBRID_DIFF_DIR/alignment
kallisto index -i AllRNASeq_Kallisto.idx $ISO_WKD_ROOT_TAMAFIL/$NAME"_sqantifiltered_tamafiltered_classification.fasta" 2> AllRNASeq_Kallisto.index.log


# run_kallisto_1sample <input_RNASEQ_rawdir> <sample> <input_ref_name_idx> <output_dir>
counter=1
for i in ${RNASEQ_SAMPLES_NAMES[@]}; do
  echo $counter
  echo $i
  run_kallisto_1sample $RNASEQ_FILTERED_DIR $i AllRNASeq_Kallisto.idx $HYBRID_DIFF_DIR/alignment
  counter=$((counter+1))
done

# generate_rnaseq_counts <input_dir>
# create one big expression matrix file from the kallisto output directory
generate_rnaseq_counts $HYBRID_DIFF_DIR