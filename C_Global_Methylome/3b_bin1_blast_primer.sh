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
#SBATCH --output=3b_bin1_blast_primer.o
#SBATCH --error=3b_bin1_blast_primer.e

# create blast library from bin1 sequences 


## print start date and time
echo Job started on:
date -u


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/C_DMP
source ${SC_ROOT}/bin1_methylation.config
source ${GENERALFUNC}/3_Transcriptome_Characterisation/blast2library.sh


cd ${BIN1_BLAST}
source activate nanopore 
Rscript $SC_ROOT/3a_bin1_blast_primer.R


##-------------------------------------------------------------------------

# generate a fasta file with bin1 specific isoforms (using ID)
Grep_Fasta ${BIN1_BLAST}/bin1_isoforms_classification.txt ${ONTFASTA} ${BIN1_BLAST}/bin1_isoforms.fasta 

# create a blast library with bin1 specific isoforms and run blast on primers
Blast_seq2seq ${BIN1_BLAST} bin1_isoforms.fasta bin1_isoforms.db ${BIN1_PRIMERS} bin1_isoforms.blast.txt build low