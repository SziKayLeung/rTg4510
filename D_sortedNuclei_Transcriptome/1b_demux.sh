#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-4429%10
#SBATCH --output=/Output/1log/1b_demux-%A_%a.o
#SBATCH --error=./Output/1log/1b_demux-%A_%a.e


# 12/02/2024: NeuN rTg4510 whole transcriptome dataset across 8 samples
# 4430 fastq files in RAW_FASTQ_2 dir

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/D_sortedNuclei_Transcriptome
source $SC_ROOT/rTg4510_snont.config
source $SC_ROOT/01_source_functions.sh


##-------------------------------------------------------------------------

raw_fastq2_files=($(ls ${RAW_FASTQ_2}/*fastq.gz))
#echo "${#raw_fastq2_files[@]}"

SamplePath=${raw_fastq2_files[${SLURM_ARRAY_TASK_ID}]}
Sample=$(basename ${SamplePath} .fastq.gz)

echo "Processing ${Sample}"

# 3) run_porechop <raw.fastq.gz> <output_dir>
mkdir -p ${WKD_ROOT}/1_demultiplex/DN/log
run_porechop ${SamplePath} ${WKD_ROOT}/1_demultiplex/DN/${Sample} > ${WKD_ROOT}/1_demultiplex/DN/log/${Sample}.log