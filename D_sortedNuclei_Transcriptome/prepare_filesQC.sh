#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=2:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --mem=200G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


# 01/2025: generate the readstats (necessary files) required to run QC from LRPipeline 

module load Miniconda2
source activate lrp 

# stats for demuxed files (raw file) 
DN=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/1b_demultiplex_merged/DN
NeuN=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/1b_demultiplex_merged/NeuN
datasets=(${DN} ${NeuN})
for dir in ${datasets[@]}; do 
	for file in ${dir}/*fastq*; do 
	cellType=$(basename $dir)
	gval=$(basename $file _merged.fastq)
	upDir=$(dirname $dir)
	echo "Processing ${file}, output: ${gval}_${cellType}_readstats.txt"
	  if [ -s ${upDir}/${gval}_${cellType}_readstats.txt ]; then
		echo "Already processed"
	  else
		echo "Running"
		seqkit stats -a ${file} > ${upDir}/${gval}_${cellType}_readstats.txt
	   fi 
	done
done


# stats for aligned files (3_minimap2)
DN=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/3_minimap/DN
NeuN=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/3_minimap/NeuN
datasets=(${DN} ${NeuN})
for dir in ${datasets[@]}; do 
	for file in ${dir}/*sorted.sam*; do 
		cellType=$(basename $dir)
		name=$(basename $file _combined_sorted.sam)
		upDir=$(dirname $dir)
		echo "Processing ${file}, output: $name"
		htsbox samview -pS $file > $upDir/${name}.paf
		awk -F'\t' '{if ($6!="*") {print $0}}' $upDir/${name}.paf > $upDir/${name}.filtered.paf
		awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $upDir/${name}.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $upDir/${name}"_mappedstats.txt"
	done
done
