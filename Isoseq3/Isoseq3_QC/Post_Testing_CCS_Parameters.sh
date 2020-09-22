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

# 05/09/2020: Post Testing_CCS_Parameters.sh with Sample K18, to run Iso-Seq3 pipeline and align with ERCC 

#************************************* DEFINE GLOBAL VARIABLES
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq
CCS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/CCS
LIMA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/LIMA
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/REFINE
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/CLUSTER
MAPPING=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/MAPPING
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/TOFU
SQANTI2_output_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/SQANTI2


#************************************* Iso-Seq3 Pipeline 
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
for i in $CCS/*ccs.bam; do 
	echo "Processing $i for Lima" 
	sample=$(basename "$i" | cut -d "_" -f 1,2,3,4 )
	run_LIMA $sample $CCS $LIMA "no_multiplex"
	run_REFINE $sample $LIMA $REFINE 
	run_CLUSTER $sample $REFINE $CLUSTER
done


#************************************* Extract Stats of CCS and Lima for downstream analysis 
module load Miniconda2 
source activate sqanti2_py3
python $FUNCTIONS/Run_Stats/CCS.py $CCS $CCS Testing
python $FUNCTIONS/Run_Stats/LIMA.py $LIMA $LIMA Testing


#************************************* Post Iso-Seq3 Pipeline
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh
for i in $CLUSTER/*hq.fasta; do 
	echo "Processing $i for Mapping" 
	sample=$(basename "$i" | cut -d "_" -f 1,2,3,4 | cut -d "." -f 1,2)
	echo $sample
	convert_fa2fq $sample $CLUSTER
	run_minimap2 $sample $CLUSTER ERCC 
	tofu $sample $CLUSTER
done 
