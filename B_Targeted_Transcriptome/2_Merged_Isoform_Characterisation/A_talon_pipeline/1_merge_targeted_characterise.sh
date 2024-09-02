#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=50:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --error=Characterise_ONT_IsoSeqMerge.e
#SBATCH --output=Characterise_ONT_IsoSeqMerge.o

# 19/01/2022: TALON, SQANTI and GFFcompare of ONT transcripts across all samples and merge with Iso-Seq

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/rTg4510_merged.config
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/01_source_function.sh
export PATH=$PATH:${LOGEN}/target_gene_annotation
export PATH=$PATH:${LOGEN}/merge_characterise_dataset

##-------------------------------------------------------------------------

#run_gff_compare <output_name> <ref_gtf> <2nd_gtf> <root_dir>
run_gff_compare $NAME $ISOSEQ_GTF $ONT_GTF $WKD_ROOT

# identify_common_transcripts <root_dir>
identify_common_transcripts $WKD_ROOT

# SQANTI annotation of merged output
run_sqanti3 $NAME"_genename" basic $WKD_ROOT

# run_cpat <input_fasta> <output_name> <output_dir>
run_cpat $WKD_ROOT/3_sqanti3/$NAME"_genename_corrected.fasta" $NAME $WKD_ROOT

# extract_best_orf <sample> <root_dir>
extract_best_orf $NAME $WKD_ROOT

# colour_by_abundance <sample> <input_gtf> <root_dir>
colour_by_abundance $NAME $WKD_ROOT/3_sqanti3/$NAME"_genename_corrected.gtf" $WKD_ROOT/2_common_transcripts/Final_Merged_Abundance.csv $WKD_ROOT

# subset_gene_reference 
subset_gene_reference $WKD_ROOT

# extract reference information
python $FICLE_ROOT/reference/extract_reference_info.py --r ${GENOME_GTF} --glist ${TGENES[@]} --split ${TGENES_REF} --exc ${TGENES_REF}/Gencode_transcript_exclusion.csv  --short_read  ${RNASEQ_COUNTS} --o ${TGENES_REF}

# run_transdecoder <name> <root_dir>
run_transdecoder $NAME $WKD_ROOT

# full_characterisation <name> <gene> <sq_filter>
# Filtered/Unfiltered refers to SQANTI filtered gtf for parsing to FICLE
Tgenes=(Rhbdf2 Clu)
for gene in ${Tgenes[@]}; do 
  full_characterisation $NAME $gene Filtered 
  full_characterisation $NAME $gene Unfiltered
done
