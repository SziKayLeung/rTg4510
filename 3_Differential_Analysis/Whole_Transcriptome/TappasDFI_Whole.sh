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
#SBATCH --reservation=research_project-mrc148213_5


#************************************* DEFINE GLOBAL VARIABLES
WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/DiffAnalysis/TAPPAS_OUTPUT/IsoSeq_Expression
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis/

module load Miniconda2/4.3.21

source activate sqanti2_py3
# R; BiocManager::install("maSigPro")

#Rscript $FUNCTIONS/TappAS_DFI.R -a0675479350 -mMASIGPRO -i$WKD/Data -d$WKD/Data/DFI -o$WKD/Data/DFI -f2.0 -tFOLD -u3 -kgroup -rLESSSTRICT -s0.05 -c$WKD/Data/Content/dfi_total_features.0675479350.tsv -g1$WKD/Data/Content/dfi_test_features.0675479350.tsv -g2$WKD/Data/Content/dfi_test_genes.0675479350.tsv -x$WKD/Data/Content/dfi_matrix.0675479350.tsv -ltNMD,CDS,5UTR_Length,polyA_Site,3UTR_Length,PAS,uORF,5UTRmotif,3UTRmotif,repeat,miRNA_Binding -lpSIGNAL,DOMAIN,ACT_SITE,INTRAMEM,TRANSMEM,PTM,BINDING,MOTIF,COMPBIAS,COILED

source deactivate 
module load R
Rscript $FUNCTIONS/TappAS_DEA.R &> TappAS_DEA.output
