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
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/MAP
# determine path directory for input data
POLISH=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/POLISH
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference

SAMPLES_NAMES=(L22 K18 S18 K17 O23)
#############################################################################################################
for map in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $map file for mapping on GMAP"
  time gmap -D $REFERENCE -d GRCm38.p4 -f samse -n 0 -t 12 -z sense_force $POLISH/"$map.polished.hq.fastq" > $map.isoforms.fastq.sam \
  2> $map.sam.log
  echo "Mapping $map succeeded"
  sort -k 3,3 -k 4,4n $map.isoforms.fastq.sam > $map.isoforms.fastq.sorted.sam
  # different sort from recommended using samtools: https://www.biostars.org/p/93559/
  samtools sort -O sam -T $map.sorted. -o $map.sort.sam $map.isoforms.fastq.sam
done

#############################################################################################################
# specific to 17th April (L22 K18 S18 K17 O23)
# 17th April: command line samtools sort as error 
#samtools sort -O sam -T L22.sorted. -o L22.sort.sam L22.isoforms.fastq.sam
#SAMPLES_NAMES=(K18 S18 K17 O23)
#for map in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $map file for mapping on GMAP"
  samtools sort -O sam -T $map.sorted. -o $map.sort.sam $map.isoforms.fastq.sam
#done
# Processing K18 file for mapping on GMAP
# Processing S18 file for mapping on GMAP
# Processing K17 file for mapping on GMAP
# Processing O23 file for mapping on GMAP

#SAMPLES_NAMES=(L22 K18 S18 K17 O23)
#for map in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $map file for mapping on GMAP"
  sort -k 3,3 -k 4,4n $map.isoforms.fastq.sam > $map.hq_isoforms.fastq.sorted.sam
#done
#Processing L22 file for mapping on GMAP
#Processing K18 file for mapping on GMAP
#Processing S18 file for mapping on GMAP
#Processing K17 file for mapping on GMAP
#Processing O23 file for mapping on GMAP


