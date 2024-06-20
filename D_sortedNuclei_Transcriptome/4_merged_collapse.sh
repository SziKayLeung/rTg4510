#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=40:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=4_merged_collapse.o
#SBATCH --error=4_merged_collapse.e


# 30minutes

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/D_sortedNuclei_Transcriptome
source $SC_ROOT/rTg4510_snont.config
source $SC_ROOT/01_source_functions.sh

LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing

export dir=$WKD_ROOT/5_cupcake
export samplename=rTg4510SCN
mkdir -p ${dir}/5_align/combined ${dir}/6_collapse ${dir}/5_align/combined_fasta ${dir}/7_sqanti3


##-------------------------------------------------------------------------

source activate nanopore
replace_filenames_with_csv.py --copy --ext=filtered.bam -i=$WKD_ROOT/5_cupcake/5_align -f=$META_ROOT/rTg4510SCN_rename.csv -d=${dir}/5_align/combined &> ${dir}/5_align/combined/rTg4510SCN_copy.log
replace_filenames_with_csv.py --copy --ext=filtered.fa -i=$WKD_ROOT/5_cupcake/5_align -f=$META_ROOT/rTg4510SCN_rename.csv -d=${dir}/5_align/combined_fasta &> ${dir}/5_align/combined_fasta/rTg4510SCN_fa_copy.log


##-------------------------------------------------------------------------

# merge alignment
echo "Collapsing..."
allfilteredmapped=($(ls ${dir}/5_align/combined/*filtered.bam)) 
for i in ${allfilteredmapped[@]}; do echo $i; done
source activate nanopore
samtools merge -f ${dir}/6_collapse/${samplename}_mapped.filtered.sorted.bam ${allfilteredmapped[@]}

# collapse
echo "Collapsing..."
echo "Output: ${dir}/6_collapse/${samplename}_collapsed.gff"
cd ${dir}/6_collapse
source activate isoseq3
isoseq3 collapse ${dir}/6_collapse/${samplename}_mapped.filtered.sorted.bam ${samplename}_collapsed.gff \
  --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
  --log-level TRACE --log-file ${samplename}_collapsed.log

# demultiplex 
source activate nanopore
adapt_cupcake_to_ont.py ${dir}/5_align/combined_fasta -o ${samplename}

demux_cupcake_collapse.py \
  ${dir}/6_collapse/${samplename}_collapsed.read_stat.txt \
  ${dir}/5_align/combined_fasta/${samplename}_sample_id.csv\
  --dataset=ont