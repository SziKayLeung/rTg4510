#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Linux_Commands_Part2.o
#SBATCH --error=Linux_Commands_Part2.e

# 11/12/2020: Run 2 samples to test output of IsoSeq and Post IsoSeq Pipeline post Tofu

#************************************* DEFINE GLOBAL VARIABLES
# File directories
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/cupcake/tofu/counting
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Testing/Linux
cd $Isoseq3_WKD
mkdir CHAIN SQANTI3
TOFU=$Isoseq3_WKD/TOFU
CHAIN=$Isoseq3_WKD/CHAIN
SQANTI3=$Isoseq3_WKD/SQANTI3

#************************************* Chain
SAMPLE_NAMES=(O18 S18)
cd $CHAIN; mkdir ${SAMPLE_NAMES[@]}

# make sure not in conda environment when renaming
for sample in ${SAMPLE_NAMES[@]}; do
	echo "Processing $sample for mm10 Chaining"
	cd $CHAIN/$sample
	cp $TOFU/{*$sample.collapsed.group.txt*,*$sample.collapsed.filtered.gff*,*$sample.collapsed.filtered.abundance.txt*,*$sample.collapsed.filtered.rep.fq*} .
	echo "Copied over files"
	ls
	rename $sample Sample *
done


# Prepare configuration file
cat << EOF > $CHAIN/Chained_Configuration.txt
SAMPLE=O18;$CHAIN/O18
SAMPLE=K18;$CHAIN/K18

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF

echo "Processing samples for chain_samples.py"
python $CUPCAKE/chain_samples.py $CHAIN/Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> Chained_Configuration.log

#************************************* SQANTI
source $FUNCTIONS/Post_Isoseq3_Function.sh
# run_sqanti3_QC <output_prefix_sample> <input_fasta> <input_abundance_file> <coverage/genome=mm10_rnqaseq/mm10/hg38_gencode/hg38_chess/ERCC> <input_kallisto_file> <input_rnaseq_dir> <output_dir>
run_sqanti3_QC chained $CHAIN/all_samples.chained.gff $CHAIN/all_samples.chained_count.txt mm10 NA NA $SQANTI3

# run_sqanti_Filter <input/output_dir> <classification.txt> <corrected.fasta> <corrected.gtf>
run_sqanti_Filter $SQANTI3 all_samples.chained_classification.txt all_samples.chained_corrected.fasta all_samples.chained_corrected.gtf
