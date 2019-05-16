#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 10th April 2019
# ran Isoseq3.sh at 11am (10th April), by 21:20, finished CCS for Sample S23_Tg4510_8months
# Error due to wrong directory for primer fasta for refine
# only rerun isoseq3 refine onwards as no error in ccs and lima
#############################################################################################################

#############################################################################################################
# Rawdata transfer
#############################################################################################################
#scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/1_C11/m54082_190403_135102.subreads.bam sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/S23
#scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/1_C11/m54082_190403_135102.subreads.bam.pbi sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/S23
#scp sl693@zeus.ex.ac.uk:/bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/1_C11/m54082_190403_135102.subreadset.xml sl693@login.isca.ex.ac.uk:/gpfs/ts0/scratch/sl693/S23
################################################################################################################

module load Anaconda2
source activate my_root 

cd /gpfs/ts0/scratch/sl693/S23/

#############################################################################################################
# Generating circular consensus sequence (ccs) from subreads
#############################################################################################################
# --minPasses=1, only need minimum 1 subread to gerenate CCS
# ccs --version
# ccs --numThreads=16 --noPolish --minPasses=1 /gpfs/ts0/scratch/sl693/S23/m54082_190403_135102.subreads.bam S23.ccs.bam

# echo ccs done

#############################################################################################################
# Isoseq3.1
#############################################################################################################
# Lima 
# removed --no pbi as this is needed for downstream polishing
# lima -version
# head /gpfs/ts0/home/sl693/reference/primer.fasta
# lima S23.ccs.bam /gpfs/ts0/home/sl693/reference/primer.fasta S23.demux.ccs.bam --isoseq --dump-clips --dump-removed --peek-guess
# echo lima done

## Isoseq3 refine from demuxed bam
isoseq3 refine --version
if isoseq3 refine S23.demux.ccs.primer_5p--primer_3p.bam /gpfs/ts0/home/sl693/reference/primer.fasta S23.flnc.bam --require-polya ; then
    echo "isoseq3 refine succeeded"
else
    echo "isoseq3 failed"
fi

# Isoseq3 cluster 
isoseq3 cluster --version
isoseq3 cluster S23.flnc.bam S23.unpolished.bam --verbose
echo cluster done

# Isoseq3 polish 
isoseq3 polish --version 
if isoseq3 polish S23.unpolished.bam m54082_190403_135102.subreadset.xml polished.bam --verbose ; then
    echo "isoseq3 polish succeeded"
else
    echo "isoseq3 polish failed"
fi
################################################################################################################
source deactivate

## command line to organise files: 16/05/2019
S23=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/S23
CCS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/CCS
LIMA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/LIMA
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/REFINE
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/CLUSTER
POLISH=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/POLISH
RAWDATA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/rawdata

cd $S23
mv ccs_report.txt S23_ccs_report.txt
mv polished.bam S23.polished.bam
mv polished.bam.pbi S23.polished.bam.pbi
mv polished.hq.fastq S23.polished.hq.fastq
mv polished.hq.fasta.gz S23.polished.hq.fasta.gz
mv polished.lq.fastq.gz S23.polished.lq.fastq.gz
mv polished.transcriptset.xml S23.polished.transcriptset.xml  
mv polished.cluster_report.csv S23.polished.cluster_report.csv
mv m54082_190403_135102.subreads* $RAWDATA
 
mv $S23/*ccs.bam* $CCS 
mv $S23/*demux* $LIMA
mv $S23/*flnc* $REFINE
mv $S23/*unpolished* $CLUSTER
mv $S23/*polished* $POLISH
mv $S23/*ccs.report* $CCS
mv m54082_190403_135102.subreads* $RAWDATA