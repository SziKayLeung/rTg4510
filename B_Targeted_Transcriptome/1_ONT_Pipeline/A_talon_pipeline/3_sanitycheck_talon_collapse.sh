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
#SBATCH --output=00_check_talon_collapse.o
#SBATCH --error=00_check_talon_collapse.e

# aim: identify the raw ont reads clustered to isoform in TALON 
# motivation: test TALON strigency for collapsing
# method:
  # convert all aligned ONT .sam to .gff3 using E.Tseng scripts
  # cat all .gff3 to one merged .gff3, representative of all reads used in TALON
  # convert .gff3 to gtf using gffread
  # run identify_raw_ont_reads.py
# nb: too memory intensive to create gff3 from one big sam file


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE/
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/
source $SC_ROOT/B_Targeted_Transcriptome/1_ONT_Pipeline/rTg4510_ont.config
source $SC_ROOT/B_Targeted_Transcriptome/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing


##-------------------------------------------------------------------------

# convert .sam to .gff3
source activate sqanti2_py3
samfiles=$(ls $WKD_ROOT/4_minimap/*/*.sam)
for i in ${samfiles[@]}; do
  echo "Processing $i" 
  sample=$(basename $i .sam)
  python ${CUPCAKE}/sequence/sam_to_gff3.py $i -s mm10 
  gffread ${sample}.gff3 -T -o ${sample}.gtf
done

# merge .gff3 to one big file
gffiles=$(ls $WKD_ROOT/4_minimap/*/*.gff3)
cat ${gffiles[@]} > $WKD_ROOT/4_minimap/all_combined_reads.gff3

# convert .gff3 to gtf using gffread
gffread all_combined_reads.gff3 -T -o all_combined_reads.gtf

# identify raw reads
# nb: also take raw_gtf from one sample (${ONT_COMBINED_RAW_GTF}) rather than combined as memory intensive
# i.e. $rTg4510_ROOT/F_ONT_Targeted/4_minimap/Batch2/BC2_combined_reads.gtf
identify_raw_talon_collapse.py ${TALON_ANNO} ${ONT_COMBINED_RAW_GTF} /
  --isoform TALONT000490198 TALONT000476505 /
  --dir=$rTg4510_ROOT/F_ONT_Targeted/6_talon/4_assess_collapse