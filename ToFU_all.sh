#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 12th April 2019
# run ToFU cupcake for all Tg4510 samples: L22, K18, O23, S18, K17 
# already Isoseq3 
#############################################################################################################
module load cDNA_Cupcake/6.9-intel-2017b-Python-2.7.14

# directory for saving output
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/ToFU
# determine path directory for input data
GMAP=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/GMAP
POLISH=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/POLISH

SAMPLES_NAMES=(L22 K18 S18 K17 O23)
#############################################################################################################
# Collapse
for sample in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $sample file for collapse"
  time collapse_isoforms_by_sam.py --input $GMAP/"$sample.polished.hq.fastq" --fq -s $GMAP/$sample.hq_isoforms.fastq.sorted.sam \
  --dun-merge-5-shorter -o $sample   
  echo "Collapse $sample successful"
done

# Create Abundance Script of full-length transcripts  
for sample in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $sample file for abundance"
  time get_abundance_post_collapse.py "$sample.collapsed" $POLISH/"$sample.cluster_report.csv"
  echo "Abundance for $sample successful"
done
 
# Remove degraded isoforms (default setting) 
for sample in "${SAMPLES_NAMES[@]}"; do 
  echo "Remove degraded isoforms for Sample $sample"
  time filter_away_subset.py "$sample.collapsed" 
  echo "Remove degraded isoforms for Sample $sample successful"
done
