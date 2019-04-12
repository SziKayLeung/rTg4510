#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 11th April 2019
# run Isoseq3 for all Tg4510 samples: L22, K18, O23, S18, K17 
# already run Isoseq3 as test on S23 sample (10th April S23_Isoseq3.sh)
#############################################################################################################
module load GMAP-GSNAP/2016-11-07-foss-2016b
gmap --version

module load SAMtools/1.7-foss-2018a
samtools --version

SAMPLES_NAMES=(L22 K18 S18 K17 O23)

cd /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/POLISH

count=0
for map in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $map file for mapping on GMAP"
  output=(${SAMPLES_NAMES[count]})
  gmap -D /gpfs/ts0/home/sl693/reference -d GRCm38.p4 -f samse -n 0 -t 12 -z sense_force "$map.polished.hq.fastq" > $output.isoforms.fastq.sam \
  2> $output.sam.log
  echo "Mapping $output succeeded"
  sort -k 3,3 -k 4,4n $output.hq_isoforms.fastq.sam > $output.hq_isoforms.fastq.sorted.sam
  # different sort from recommended using samtools: https://www.biostars.org/p/93559/
  samtools sort -O sam -T $output.sorted. -o $output.sort.sam $output.hq_isoforms.fastq.sam
  count=$((count+1)) 
done