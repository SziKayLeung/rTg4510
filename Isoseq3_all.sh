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
module load Anaconda2
source activate my_root 

# Listing versions 
ccs --version
lima --version 
isoseq3 --version
#############################################################################################################
module load Anaconda2
source activate my_root 

head samples.txt

SAMPLES_NAMES=(L22 K18 S18 K17 O23)

cd /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/Parameters
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

#############################################################################################################
# Isoseq3.1 
#############################################################################################################

# Lima 
# removed --no pbi as this is needed for downstream polishing
for lima in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $lima file for demultiplexing"
  time lima "$lima.ccs.bam" $FASTA $lima.demux.ccs.bam --isoseq --dump-clips --dump-removed --peek-guess
  echo "lima $lima succeeded"
done

## Isoseq3 refine from demuxed bam
for refine in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $refine file for refine"
  time isoseq3 refine "$refine.demux.ccs.primer_5p--primer_3p.bam" $FASTA $refine.flnc.bam --require-polya
  echo "refine $refine succeeded"
done

# Isoseq3 cluster 
for cluster in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $cluster file for cluster"
  time isoseq3 cluster "$cluster.flnc.bam" $cluster.unpolished.bam --verbose
  echo "cluster $cluster succeeded"
done

ls *unpolished.bam
# Isoseq3 polish 
count=0
for polish in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $polish file..."
  time isoseq3 polish "$polish.unpolished.bam" ${SUB_FILES[count]} $polish.polished.bam --verbose    
  echo "polish $polish succeeded"
  gunzip $polish.polished.hq.fastq
  echo "unzipped $polish.polished.hq.fastq successful"  
  count=$((count+1))
done
#############################################################################################################
source deactivate
mv *ccs.bam* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/CCS
mv *demux* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/LIMA
mv *flnc* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/REFINE
mv *unpolished* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/CLUSTER
mv *polished* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/POLISH
mv *ccs.report* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/CCS

echo done

