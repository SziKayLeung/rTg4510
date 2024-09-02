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
#SBATCH --array=0-7%4 # 8 barcodes per flow cell
#SBATCH --output=Output/3log/3_pbalign_filter_DN-%A_%a.o
#SBATCH --error=Output/3log/3_pbalign_filter_DN-%A_%a.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/D_sortedNuclei_Transcriptome
source $SC_ROOT/rTg4510_snont.config
source $SC_ROOT/01_source_functions.sh

DNfa=($(find "${WKD_ROOT}/4_tclean/DN" -type f -name "*clean.fa")) 
fasta=${DNfa[${SLURM_ARRAY_TASK_ID}]}
sample=$(basename ${fasta} | cut -d "_" -f 1)


##-------------------------------------------------------------------------

# align
echo "Aligning ${sample}: ${fasta} ..."
echo "Output: ${WKD_ROOT}/5_cupcake/5_align/${sample}_DN_mapped.bam"
source activate isoseq3
cd ${WKD_ROOT}/5_cupcake/5_align
pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${fasta} ${sample}_DN_mapped.bam --log-level TRACE --log-file ${sample}_mapped.log

# filter_alignment <input_name> <input_mapped_dir>
# output = ${sample}_mapped.filtered.bam, ${sample}_mapped.filtered.sorted.bam
filter_alignment ${sample}_DN_mapped ${WKD_ROOT}/5_cupcake/5_align