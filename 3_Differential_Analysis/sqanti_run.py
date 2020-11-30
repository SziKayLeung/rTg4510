#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:30:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
SQANTI3_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/SQANTI3
CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence

module load Miniconda2
source activate sqanti2_py3 
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual_Samples/CHAIN/SQANTI3
export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf ../all_samples.chained.gff $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --fl_count ../all_samples.chained_count.txt
