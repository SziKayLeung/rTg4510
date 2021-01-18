#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Individual_Part2.o
#SBATCH --error=Individual_Part2.e

# 13/12/2020: Repeat Chaining with updated cupcake v17.0. (reinstalled on SQANTI2_py3)

FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
SQANTI3=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/SQANTI
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/TOFU
CHAIN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/CHAIN
SUBCHAIN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/SUBCHAIN

cp -R $TOFU/* $CHAIN/
cp -R $TOFU/* $SUBCHAIN/

cat << EOF > $CHAIN/Chained_Configuration.txt
SAMPLE=M21;$CHAIN/M21
SAMPLE=K23;$CHAIN/K23
SAMPLE=O18;$CHAIN/O18
SAMPLE=K18;$CHAIN/K18
SAMPLE=S18;$CHAIN/S18
SAMPLE=L22;$CHAIN/L22
SAMPLE=Q20;$CHAIN/Q20
SAMPLE=K24;$CHAIN/K24
SAMPLE=Q21;$CHAIN/Q21
SAMPLE=K17;$CHAIN/K17
SAMPLE=O23;$CHAIN/O23
SAMPLE=S23;$CHAIN/S23

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF

cat << EOF > $SUBCHAIN/Chained_Configuration.txt
SAMPLE=M21;$CHAIN/M21
SAMPLE=K23;$CHAIN/K23
SAMPLE=S23;$CHAIN/S23

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF


module load Miniconda2
source activate cupcake

cd $SUBCHAIN
chain_samples.py $SUBCHAIN/Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> Chained_Configuration.log

cd $CHAIN
echo "Processing samples for chain_samples.py"
chain_samples.py $CHAIN/Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> Chained_Configuration.log


#************************************* SQANTI
source $FUNCTIONS/Post_Isoseq3_Function.sh
# run_sqanti3_QC <output_prefix_sample> <input_fasta> <input_abundance_file> <coverage/genome=mm10_rnqaseq/mm10/hg38_gencode/hg38_chess/ERCC> <input_kallisto_file> <input_rnaseq_dir> <output_dir>
# run_sqanti_Filter <input/output_dir> <classification.txt> <corrected.fasta> <corrected.gtf>
#run_sqanti3_QC chained $CHAIN/all_samples.chained.gff $CHAIN/all_samples.chained_count.txt mm10 NA NA $SQANTI3
#run_sqanti_Filter $SQANTI3 all_samples.chained_classification.txt all_samples.chained_corrected.fasta all_samples.chained_corrected.gtf

#for i in ${SAMPLES_NAMES[@]}; do
#	run_sqanti3_QC chained $TOFU/$i".collapsed.filtered.gff" $TOFU/$i."collapsed.filtered.abundance.txt" mm10 NA NA $SQANTI3
	#run_sqanti_Filter $SQANTI3 $i".collapsed.filtered_classification.txt" $i".collapsed.filtered_corrected.fasta" $i".collapsed.filtered_corrected.gtf"
#done
