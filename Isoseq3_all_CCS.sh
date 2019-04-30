#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 16th April 2019
# rerun CCS for all Tg4510 samples as previously missed out on ccs report: L22, K18, O23, S18, K17 
#############################################################################################################
module load Anaconda2
source activate my_root 

# Listing versions 
ccs --version
lima --version 
isoseq3 --version
#############################################################################################################
module load Anaconda2
source activate my_root 

SAMPLES_NAMES=(L22 K18 S18 K17 O23)

cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/Parameters
# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
head raw.txt 
BAM_FILES=(`cat "raw.txt" | egrep -v "^\s*(#|$)"`)

head sub.txt
SUB_FILES=(`cat "sub.txt" | egrep -v "^\s*(#|$)"`)

FASTA=/gpfs/ts0/home/sl693/reference/primer.fasta
head $FASTA

# change directory for saving 
cd /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3
#############################################################################################################
# Generating circular consensus sequence (ccs) from subreads
#############################################################################################################

count=0
for f in "${BAM_FILES[@]}"; do
  echo "Processing $f file..."
  output=(${SAMPLES_NAMES[count]})
  time ccs --numThreads=16 --noPolish --minPasses=1 $f $output.ccs.bam
  mv ccs_report.txt $output.ccs.report.txt
  count=$((count+1)) 
done
ls *ccs.bam

mv *ccs.report* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/CCS