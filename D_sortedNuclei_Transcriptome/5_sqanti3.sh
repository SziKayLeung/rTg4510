#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=5_sqanti3.o
#SBATCH --error=5_sqanti3.e


# Batch2: realign with pbmm2 align and filter alignment 

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

# sqanti3
echo "Running SQANTI3..."
source activate sqanti2_py3
cd ${dir}/7_sqanti3
python $SQANTI3_DIR/sqanti3_qc.py ${dir}/6_collapse/${samplename}_collapsed.gff \
$GENOME_GTF $GENOME_FASTA -t 30 --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA \
--genename --isoAnnotLite --skipORF --report skip &> ${samplename}_sqanti_qc.log

# sqanti3 filter 
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_collapsed_classification.txt \
--faa=${samplename}_collapsed_corrected.fasta \
--gtf=${samplename}_collapsed_corrected.gtf \
-j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log