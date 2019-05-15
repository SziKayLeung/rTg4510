#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 11th April 2019
# 11-14th April: run Isoseq3 for all Tg4510 samples: L22, K18, O23, S18, K17 
# 15th May: run Isoseq3 for additional WT samples for paper: C21
# not able to access raw data for sample C20, Q21
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

# File directories 
PARAMETERS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/Parameters
CCS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/CCS
LIMA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/LIMA
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/REFINE
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/CLUSTER
POLISH=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/POLISH

SAMPLES_NAMES=(L22 K18 S18 K17 O23 C21)
# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
cd $PARAMETERS
cat raw.txt 
BAM_FILES=(`cat "raw.txt" | egrep -v "^\s*(#|$)"`)

cat sub.txt
SUB_FILES=(`cat "sub.txt" | egrep -v "^\s*(#|$)"`)

cat primer.fasta
FASTA=$PARAMETERS/primer.fasta

#############################################################################################################
# Generating circular consensus sequence (ccs) from subreads
#############################################################################################################

count=0
cd $CCS
for f in "${BAM_FILES[@]}"; do
  echo "Processing $f"
  output=(${SAMPLES_NAMES[count]})
  echo "Checking if $output.ccs.bam exists"
  if [ -f $output.ccs.bam ]; then 
    echo "$output.ccs.bam file already exists; CCS no need to be processed on Sample $output"
  else 
    echo "$output file does not exist"
    echo "Processing CCS for sample $output"
    time ccs --numThreads=16 --noPolish --minPasses=1 $f $output.ccs.bam
    echo "CCS for Sample $output successful"
    mv ccs_report.txt $output.ccs.report.txt
  fi
  count=$((count+1)) 
done 

#############################################################################################################
# Isoseq3.1.2 
#############################################################################################################

# Lima 
# removed --no pbi as this is needed for downstream polishing
cd $LIMA
for lima in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $lima file for demultiplexing"
  if [ -f $lima.demux.ccs.json ]; then 
    echo "$lima.demux.ccs.json file already exists; LIMA no need to be processed on Sample $lima"
  else 
    echo "$lima.demux.ccs.bam file does not exist"
    time lima "$lima.ccs.bam" $FASTA $lima.demux.ccs.bam --isoseq --dump-clips --dump-removed --peek-guess
    echo "lima $lima successful"
  fi
done

## Isoseq3 refine from demuxed bam
cd $REFINE
for refine in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $refine file for refine"
  if [ -f $refine.flnc.bam ]; then
    echo "$refine.flnc bam file already exists; Refine no need to be processed on Sample $refine"
  else
    echo "$refine.flnc bam file does not exist"
    time isoseq3 refine "$refine.demux.ccs.primer_5p--primer_3p.bam" $FASTA $refine.flnc.bam --require-polya
    echo "refine $refine successful"
  fi
done

# Isoseq3 cluster
cd $CLUSTER
for cluster in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $cluster file for cluster"
  if [ -f $cluster.unpolished.bam ]; then
    echo "$cluster.unpolished.bam file already exists; Cluster no need to be processed on Sample $cluster"
  else 
    echo "$cluster.unpolished.bam file does not exist"
    #time isoseq3 cluster "$cluster.flnc.bam" $cluster.unpolished.bam --verbose
    echo "cluster $cluster successful"
  fi
done

ls *unpolished.bam

# Isoseq3 polish 
cd $POLISH
count=0
for polish in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $polish file..."
  if [ -f $polish.polished.hq.fasta ]; then
    echo "$polish.polished.hq.fasta file already exists; Polish no need to be processed on Sample $polish"
  else 
    echo "$polish.polished.hq.fasta file does not exist"
    #time isoseq3 polish "$polish.unpolished.bam" ${SUB_FILES[count]} $polish.polished.bam --verbose    
    echo "polish $polish successful"
    #gunzip $polish.polished.hq.fastq
    echo "unzipped $polish.polished.hq.fastq successful" 
  fi 
  count=$((count+1))
done
#############################################################################################################
source deactivate
# first processing of bash loop filing; 
#mv *ccs.bam* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/CCS
#mv *demux* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/LIMA
#mv *flnc* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/REFINE
#mv *unpolished* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/CLUSTER
#mv *polished* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/POLISH
#mv *ccs.report* /gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/CCS

echo IsoSeq3 done

