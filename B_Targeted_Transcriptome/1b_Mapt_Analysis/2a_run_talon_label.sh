#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=100:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-8 # 9 samples
#SBATCH --output=TalonLabel-%A_%a.o
#SBATCH --error=TalonLabel-%A_%a.e

#************************************* DEFINE GLOBAL VARIABLES
WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON_Human
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019 # directory of genome files

SAMPLES=(BC1 BC2 BC3 BC4 BC5 BC6 BC7 BC8 BC9) # 9 samples for 2 batches == 18 samples
Sample=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}


module load Miniconda2
source activate sqanti2_py3

echo "Label ${Sample} for TALON"
talon_label_reads --f $WKD/Mapped_reads/Batch2/${Sample}"_clean.sam" --g $REFERENCE/hg38.fa --t 1 --ar 20 --tmpDir=$WKD/Labelled/Batch2/${Sample}_label_reads --deleteTmp --o $WKD/Labelled/Batch2/${Sample}
talon_label_reads --f $WKD/Mapped_reads/Batch3/${Sample}"_clean.sam" --g $REFERENCE/hg38.fa --t 1 --ar 20 --tmpDir=$WKD/Labelled/Batch3/${Sample}_label_reads --deleteTmp --o $WKD/Labelled/Batch3/${Sample}

