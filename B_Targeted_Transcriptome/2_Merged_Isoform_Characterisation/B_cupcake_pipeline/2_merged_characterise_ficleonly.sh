#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


# 22/01/2022: Rerun FICLE due to gfps isca down

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2
source activate ficle
FICLE_ROOT=/lustre/projects/Research_Project-MRC148213/sl693/scripts/FICLE
export PATH=$PATH:${FICLE_ROOT}


##-------------------------------------------------------------------------

# directory paths
TGENES=(Apoe Clu App Snca Ptk2b Bin1 Fus Vgf Picalm Mapt Trem2 Tardbp Sorl1 Abca7 Fyn Abca1 Cd33 Ank1 Rhbdf2 Trpa1)
refDir=/lustre/projects/Research_Project-MRC148213/sl693/reference/annotation/
inputDir=/lustre/projects/Research_Project-MRC148213/sl693/rTg4510/G_Merged_Targeted/2_sqanti3/
outputDir=/lustre/projects/Research_Project-MRC148213/sl693/rTg4510/G_Merged_Targeted/4_characterise/
  
# run ficle
for g in ${TGENES[@]}; do
  echo $g
  ficle.py --gene=$g \
    --reference=${refDir}/gencode.M22.annotation.20Targets.gtf \
    --input_bed=${inputDir}/all_iso_ont_collapsed.filtered_counts_filtered_sorted.bed12 \
    --input_gtf=${inputDir}/all_iso_ont_collapsed.filtered_counts_filtered.gtf  \
    --input_class=${inputDir}/all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt  \
    --output_dir=${outputDir} &> ${outputDir}/Log/$g"_characterise.log"
done
