#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Individual_Part3.o
#SBATCH --error=Individual_Part3.e

# 14/12/2020: TAMA merge output from collapsed.gtf

#************************************* DEFINE GLOBAL VARIABLES
Individual=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/
TAMA_prepare_rscript=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/TAMA/TAMA_Merge_Prepare.R
TAMA_Filelist=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Whole_Transcriptome/Individual_Tama_filelist.txt
TAMA_software=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama

module load Miniconda2
# Modify_genomegtf_TAMAinput <input/output_name> <input/output_dir>
Modify_gtf_TAMAinput(){

	source activate nanopore
	###################################
	# Convert genome gtf to bed12
	###################################
	cd $2
	gtfToGenePred $1.collapsed.gff $1.genepred
	genePredToBed $1.genepred $1.bed12

	source deactivate
}

#SAMPLES_NAMES=(O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23)
#for i in ${SAMPLES_NAMES[@]}; do Modify_gtf_TAMAinput $i $Individual/TOFU; done
# Rscript script.R <name>.bed12 <input_dir>
#for i in ${SAMPLES_NAMES[@]}; do Rscript $TAMA_prepare_rscript $i $Individual/TOFU; done

# prepare TAMA config file
#DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual_Samples_OLD/SQANTI
#for i in O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23; do
#  echo "$DIR/$i"_mod.bed"12:no_cap:1,1,1:$i"|sed s/":"/"\t"/g>>$TAMA_Filelist
#done

source activate sqanti2
cd $Individual/TAMA
python $TAMA_software/tama_merge.py -f $TAMA_Filelist -a 50 -z 50 -m 20 -p Individual &> Tama_merge.log
source deactivate
