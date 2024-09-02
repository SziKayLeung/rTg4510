#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=1a11_demux.o
#SBATCH --error=1a11_demux.e

# 12/02/2024: NeuN rTg4510 whole transcriptome dataset across 8 samples
# 7775 fastq files in RAW_FASTQ_1 dir
# missing samples after processing the first 10 batches

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/D_sortedNuclei_Transcriptome
source $SC_ROOT/rTg4510_snont.config
source $SC_ROOT/01_source_functions.sh


##-------------------------------------------------------------------------

missingSamplesFile=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/1_demultiplex/NeuN/missingSamples.txt
missingSamples=($(awk -F "\"*,\"*" '{print $1}' ${missingSamplesFile}))

for Sample in ${missingSamples[@]}; do 
  SamplePath=${RAW_FASTQ_1}/${Sample}.fastq.gz
  echo "Processing ${Sample}"
  #ls ${SamplePath}
  
  # 3) run_porechop <raw.fastq.gz> <output_dir>
  #mkdir -p ${WKD_ROOT}/1_demultiplex/NeuN/log
  run_porechop ${SamplePath} ${WKD_ROOT}/1_demultiplex/NeuN/${Sample} > ${WKD_ROOT}/1_demultiplex/NeuN/log/${Sample}.log
done