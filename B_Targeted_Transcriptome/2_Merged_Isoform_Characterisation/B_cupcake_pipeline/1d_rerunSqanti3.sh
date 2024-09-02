#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=1d_rerunSqanti3.o
#SBATCH --error=1d_rerunSqanti3.e


##-------------------------------------------------------------------------

# 05/02/2024: Align RNA-Seq to all_iso_ont collapsed data and re-run SQANTI3

module load Miniconda2/4.3.21
FICLE_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/FICLE/
LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510
source $SC_ROOT/B_Targeted_Transcriptome/2_Merged_Isoform_Characterisation/rTg4510_merged.config

NAME="all_iso_ont"
OUTPUTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/B_cupcake_pipeline/3_sqanti3/reRunValidation
RNASEQ_MAPPED_DIR=${RNASEQ_MAPPED_DIR}/All

##-------------------------------------------------------------------------

echo "Running SQANTI3..."
source activate sqanti2_py3
cd ${OUTPUTDIR}
python $SQANTI3_DIR/sqanti3_qc.py $CUPMERGE_DIR/2_collapse/${NAME}_collapsed.gff \
$GENOME_GTF $GENOME_FASTA -t 30 --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA -c "${RNASEQ_MAPPED_DIR}/*SJ.out.bed" \
--genename --skipORF --report skip &> ${NAME}_sqanti_qc.log

# sqanti3 filter 
python $SQANTI3_DIR/sqanti3_filter.py rules ${NAME}_collapsed_classification.txt \
--faa=${NAME}_collapsed_corrected.fasta \
--gtf=${NAME}_collapsed_corrected.gtf \
-j=${filteringJson} --skip_report &> ${NAME}_sqanti_filter.log