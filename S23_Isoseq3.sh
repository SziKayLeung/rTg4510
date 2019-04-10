#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 10th April 2019
# CCS for Sample S23_Tg4510_8months
# Isoseq3.1
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
ccs --version
ccs --numThreads=16 --noPolish --minPasses=1 /gpfs/ts0/scratch/sl693/S23/m54082_190403_135102.subreads.bam S23.ccs.bam

echo ccs done

#############################################################################################################
# Isoseq3.1
#############################################################################################################
# Lima 
# removed --no pbi as this is needed for downstream polishing
lima -version
head /gpfs/ts0/home/sl693/reference/primer.fasta
lima S23.ccs.bam /gpfs/ts0/home/sl693/reference/primer.fasta S23.demux.ccs.bam --isoseq --dump-clips --dump-removed --peek-guess
echo lima done

## Isoseq3 refine from demuxed bam
isoseq3 refine --version
isoseq3 refine S23.demux.ccs.primer_5p--primer_3p.bam /mnt/data1/Szi/reference/primer.fasta S23.flnc.bam --require-polya 
echo refine done

# Isoseq3 cluster 
isoseq3 cluster --version
isoseq3 cluster S23.flnc.bam S23.unpolished.bam --verbose -j 32
echo cluster done

# Isoseq3 polish 
isoseq3 polish --version 
isoseq3 polish S23.unpolished.bam m54082_190403_135102.subreadset.xml polished.bam --verbose -j 32
echo polish done
################################################################################################################
source deactivate
