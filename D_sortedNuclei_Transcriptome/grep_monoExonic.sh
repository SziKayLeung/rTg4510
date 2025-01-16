#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=0:30:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.k.leung@exeter.ac.uk # email address


# 01/2025: to find the number of mono-exonic reads in the aligned bam files

##-------------------------------------------------------------------------

collapseDir=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/5_cupcake/6_collapse
sqantiDir=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/5_cupcake/7_sqanti3
alignedDir=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/5_cupcake/5_align

##-------------------------------------------------------------------------

module load Miniconda2
source activate nanopore 

for bam in ${alignedDir}/*_mapped.bam; do 
	echo "For sample: ${bam}"
	total_reads=$(samtools view -c -F 4 ${bam})
	monoexonic_reads=$(samtools view -F 4 ${bam} | awk '$6 !~ /N/ {monoexonic++} END {print monoexonic}')
	percentage=$(echo "scale=2; $monoexonic_reads / $total_reads * 100" | bc)
	echo "Percentage of monoexonic reads: $percentage%"
done
