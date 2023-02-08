#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --error=4_collapse_further.e
#SBATCH --output=4_collapse_further.o


# Jan 2023: run Cupcake collapse (v8.6) on merged Iso-Seq+ONT dataset
# Premise:
  # noticed many redundant isoforms from TALON that can be further collapsed
# Method:
  # run minimap2 again for input into cupcake collapse
  # run cupcake collapse with different parameters for --5_diff and 3_diff of exon ends
  # identify the raw reads that were collapsed under the parameters (for visualisation)
  # generate the abundance file 


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/rTg4510_merged.config
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/01_source_function.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing


##-------------------------------------------------------------------------

source activate sqanti2_py3

# run_minimap <output_name> <output_dir> <input_fasta>
run_minimap IsoSeqONT $WKD_ROOT/4_characterise/CollapseMore $WKD_ROOT/3_sqanti3/IsoSeqONT.collapsed_corrected.fasta

# optimise_collapse <output_name> <input_fasta> <input_map_dir> <max_5_diff> <max_3_diff> 
max_3_diff=(10 5 20 50 100)
max_5_diff=(10 5 20 50 1000)
for i in {2..4}; do 
  echo "******************* max3diff = ${max_3_diff[$i]} and max5diff = ${max_5_diff[$i]}"
  optimise_collapse ${NAME} $WKD_ROOT/3_sqanti3/IsoSeqONT.collapsed_corrected.fasta $WKD_ROOT/4_characterise/CollapseMore/map ${max_5_diff[$i]} ${max_3_diff[$i]} 
done

# identify raw reads for collapsing to visualise on UCSC genome browser
dir=$MERGE_DIR/4_characterise/CollapseMore
inGtf=$WKD_ROOT/3_sqanti3/IsoSeqONT.collapsed_corrected.gtf
identify_raw_cupcake_collapse.py $dir/50_5diff_50_3diff/IsoSeqONT.collapsed.group.txt $inGtf --isoform TALONT000440029 --rep
identify_raw_cupcake_collapse.py $dir/1000_5diff_100_3diff/IsoSeqONT.collapsed.group.txt $inGtf --isoform TALONT000440029 --rep
identify_raw_cupcake_collapse.py $dir/1000_5diff_100_3diff/IsoSeqONT.collapsed.group.txt $inGtf --isoform TALONT000440029 --rep  

# generate abundance file
demux_cupcake_collapse_v8.6.py ${MERGE_COUNTS} ${dir}/20_5diff_20_3diff/IsoSeqONT.collapsed.group.txt ${MERGE_SAMPLES_EXP} 