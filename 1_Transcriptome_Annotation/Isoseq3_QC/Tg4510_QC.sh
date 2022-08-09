#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


# 03/09/2020: SequelTools and Lengths QC for Tg4510 mouse 
 
#************************************* DEFINE GLOBAL VARIABLES
POLISHED=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/CLUSTER
LENGTHS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/Isoseq3_WKD/Lengths
SEQUELTOOLS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/SequelTools/Scripts
RAWDATA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/rawdata

#************************************* Polished Lengths 
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq/Run_Stats/FastqLength_Functions.sh
SAMPLES=(B21 E18 K23 M21 Q20 S23 C20 K17 K24 O18 Q21 C21 K18 L22 O23 S18)

#for i in ${SAMPLES[@]}; do 
	#polished_hq_length $i $POLISHED $LENGTHS
#done


#************************************* Sequel QC Tools 
module load Miniconda2 
source activate sqanti2_py3

cd $SEQUELTOOLS
export PATH=$PATH:"$(pwd)"

cd $RAWDATA/QC
find $RAWDATA -name "*subreads.bam" > subreads.txt
bash SequelTools.sh -t Q -u subreads.txt -p a -k
