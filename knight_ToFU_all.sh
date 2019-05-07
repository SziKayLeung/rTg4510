#!/bin/sh

##############################################################################################################
# Date: 7th May 2019
# run ToFU cupcake for all Tg4510 samples: L22, K18, O23, S18, K17 
# already Isoseq3 
#############################################################################################################

# directory for saving output
cd /mnt/data1/Szi/IsoSeq_WholeTranscriptome/ToFU
# determine path directory for input data
GMAP=/mnt/data1/Szi/IsoSeq_WholeTranscriptome/Mapped
POLISH=/mnt/data1/Szi/IsoSeq_WholeTranscriptome/Isoseq/POLISH
COLLAPSE=/home/sLeung/cDNA_Cupcake/cupcake/tofu

SAMPLES_NAMES=(L22 K18 S18 K17 O23)
source activate anaCogent
#############################################################################################################
# Collapse
for sample in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $sample file for collapse"
  if time python $COLLAPSE/collapse_isoforms_by_sam.py --input $POLISH/"$sample.polished.hq.fastq" --fq -s $GMAP/$sample.hq_isoforms.fastq.sorted.sam \
  --dun-merge-5-shorter -o $sample ; then
      echo "Collapse $sample successful"
  else
      echo "Collapse $sample failed"
   fi      
done

# Create Abundance Script of full-length transcripts  
for sample in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $sample file for abundance"
  if time get_abundance_post_collapse.py "$sample.collapsed" $POLISH/"$sample.polished.cluster_report.csv" ; then 
      echo "Abundance for $sample successful"
  else 
      echo "Abundance for $sample failed"
  fi
done
 

# Remove degraded isoforms (default setting) 
for sample in "${SAMPLES_NAMES[@]}"; do 
  echo "Remove degraded isoforms for Sample $sample"
  if time filter_away_subset.py "$sample.collapsed" ; then
      echo "Remove degraded isoforms for Sample $sample successful"
  else 
      echo "Remove degraded isoforms for Sample $sample failed"
  fi 
done
conda deactivate