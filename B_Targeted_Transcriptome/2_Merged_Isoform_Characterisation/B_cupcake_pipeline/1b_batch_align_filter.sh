#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-18%9 #19 samples (9 ONT batch2, 9 ONT batch 3, 1 merged Iso-Seq)
#SBATCH --output=../../../bash_output/1b_merge_targeted-%A_%a.o
#SBATCH --error=../../../bash_output/1b_merge_targeted-%A_%a.e

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/rTg4510_merged.config
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/01_source_function.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous

##-------------------------------------------------------------------------

export dir=${CUPMERGE_DIR}

all_ont_fa=($(ls ${ONT_TCLEAN_DIR}/AllBatch2Batch3/*.fa))  

cp ${ISO_MERGED_CLUSTER_DIR}/AllMouseTargeted.clustered.hq.fasta ${dir}/1_merge_collapse/AllIsoTargeted_clustered.hq.fasta 
iso_fa=${dir}/1_merge_collapse/AllIsoTargeted_clustered.hq.fasta 

all_iso_ont_fa=(${all_ont_fa[@]} ${iso_fa})

fasta=${all_iso_ont_fa[${SLURM_ARRAY_TASK_ID}]}
sample=$(basename ${fasta} | cut -d "_" -f 1)

##-------------------------------------------------------------------------

# align
echo "Aligning ${sample}: ${fasta} ..."
echo "Output: ${dir}/1_align/${sample}_mapped.bam"
source activate isoseq3
cd ${dir}/1_align
pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} ${fasta} ${sample}_mapped.bam --log-level TRACE --log-file ${sample}_mapped.log

# filter_alignment <input_name> <input_mapped_dir>
# output = ${sample}_mapped.filtered.bam, ${sample}_mapped.filtered.sorted.bam
filter_alignment ${sample}_mapped ${dir}/1_align