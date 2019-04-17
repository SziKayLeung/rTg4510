#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=10:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 12th April 2019
# run GMAP for all Tg4510 samples: L22, K18, O23, S18, K17 
# already Isoseq3
#############################################################################################################
module load GMAP-GSNAP/2018-07-04-foss-2018b
gmap --version
module load SAMtools/1.7-foss-2018a
samtools --version
#############################################################################################################
# directory for saving output
cd /gpfs/ts0/scratch/sl693/WholeTranscriptome/MAP
# determine path directory for input data
POLISH=/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/POLISH

SAMPLES_NAMES=(L22 K18 S18 K17 O23)
#############################################################################################################
for map in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $map file for mapping on GMAP"
  time gmap -D /gpfs/ts0/home/sl693/reference -d GRCm38.p4 -f samse -n 0 -t 12 -z sense_force $POLISH/"$map.polished.hq.fastq" > $map.isoforms.fastq.sam \
  2> $map.sam.log
  echo "Mapping $map succeeded"
  sort -k 3,3 -k 4,4n $map.hq_isoforms.fastq.sam > $map.hq_isoforms.fastq.sorted.sam
  # different sort from recommended using samtools: https://www.biostars.org/p/93559/
  samtools sort -O sam -T $map.sorted. -o $map.sort.sam $map.hq_isoforms.fastq.sam
done