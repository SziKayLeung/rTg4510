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
#SBATCH --array=0-2 # 3 samples
#SBATCH --output=Post_IsoSeq3b-%A_%a.o
#SBATCH --error=Post_IsoSeq3b-%A_%a.e
# 22/10/2020: Merging All, TG and WT samples after Isoseq3_all.sh

#************************************* DEFINE GLOBAL VARIABLES
# File directories 
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/Isoseq3_WKD
REFINE=$Isoseq3_WKD/REFINE
All_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
TG_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/TG_Merged
WT_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/WT_Merged

#************************************* TO RUN FUNCTIONS ON WORKING SCRIPT
source $FUNCTIONS/Isoseq3.2.2_Functions.sh

#************************************* All samples, WT, TG merged at refine  
# merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
merging_at_refine $REFINE $All_Isoseq3_WKD/CLUSTER All_Merged O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23
merging_at_refine $REFINE $TG_Isoseq3_WKD/CLUSTER TG_Merged O18 K18 S18 L22 Q20 K24
merging_at_refine $REFINE $WT_Isoseq3_WKD/CLUSTER WT_Merged Q21 K17 M21 O23 S23 K23

#************************************* DEFINE GLOBAL VARIABLES
CHAIN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/CHAINED_ANALYSIS
CHAIN_INDIVIDUAL=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/CHAIN

RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
KALLISTO=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/Kallisto
All_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
TG_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/TG_Merged
WT_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/WT_Merged
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh

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
#SBATCH --array=0-2 # 3 samples
#SBATCH --output=Post_IsoSeq3b-%A_%a.o
#SBATCH --error=Post_IsoSeq3b-%A_%a.e

# 23/10/2020: Minimap2, TOFU Merged, WT vs TG samples 

#************************************* DEFINE GLOBAL VARIABLES
All_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
TG_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/TG_Merged
WT_Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/WT_Merged

source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh

#cd $All_Isoseq3_WKD; mkdir MAPPING TOFU; cd MAPPING; mkdir ERCC mm10; cd .. ; cd TOFU; mkdir ERCC mm10 
#cd $WT_Isoseq3_WKD; mkdir MAPPING TOFU; cd MAPPING; mkdir ERCC mm10; cd .. ; cd TOFU; mkdir ERCC mm10 
#cd $TG_Isoseq3_WKD; mkdir MAPPING TOFU; cd MAPPING; mkdir ERCC mm10; cd .. ; cd TOFU; mkdir ERCC mm10 

### setting array names for samples
SAMPLES_BATCHES=(All_Merged TG_Merged WT_Merged) 
SAMPLES_BATCHES_CLUSTER_DIR=($All_Isoseq3_WKD/CLUSTER $TG_Isoseq3_WKD/CLUSTER $WT_Isoseq3_WKD/CLUSTER)
SAMPLES_BATCHES_MM10_MAPPING_DIR=($All_Isoseq3_WKD/MAPPING/mm10 $TG_Isoseq3_WKD/MAPPING/mm10 $WT_Isoseq3_WKD/MAPPING/mm10)
SAMPLES_BATCHES_ERCC_MAPPING_DIR=($All_Isoseq3_WKD/MAPPING/ERCC $TG_Isoseq3_WKD/MAPPING/ERCC $WT_Isoseq3_WKD/MAPPING/ERCC)
SAMPLES_BATCHES_MM10_TOFU_DIR=($All_Isoseq3_WKD/TOFU/mm10 $TG_Isoseq3_WKD/TOFU/mm10 $WT_Isoseq3_WKD/TOFU/mm10)
SAMPLES_BATCHES_ERCC_TOFU_DIR=($All_Isoseq3_WKD/TOFU/ERCC $TG_Isoseq3_WKD/TOFU/ERCC $WT_Isoseq3_WKD/TOFU/ERCC)

SAMPLES=${SAMPLES_BATCHES[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_CLUSTER_DIR=${SAMPLES_BATCHES_CLUSTER_DIR[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_MM10_MAPPING_DIR=${SAMPLES_BATCHES_MM10_MAPPING_DIR[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_ERCC_MAPPING_DIR=${SAMPLES_BATCHES_ERCC_MAPPING_DIR[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_MM10_TOFU_DIR=${SAMPLES_BATCHES_MM10_TOFU_DIR[${SLURM_ARRAY_TASK_ID}]}
SAMPLES_ERCC_TOFU_DIR=${SAMPLES_BATCHES_ERCC_TOFU_DIR[${SLURM_ARRAY_TASK_ID}]}

#************************************* All samples, WT, TG Minimap 

# convert_fa2fq <file_name> <input_dir>
convert_fa2fq All_Merged.clustered.hq.fasta $All_Isoseq3_WKD/CLUSTER
convert_fa2fq WT_Merged.clustered.hq.fasta $WT_Isoseq3_WKD/CLUSTER
convert_fa2fq TG_Merged.clustered.hq.fasta $TG_Isoseq3_WKD/CLUSTER

# run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
#run_minimap2 ${SAMPLES} ${SAMPLES_CLUSTER_DIR} mm10 ${SAMPLES_MM10_MAPPING_DIR} 
#run_minimap2 ${SAMPLES} ${SAMPLES_CLUSTER_DIR} ERCC ${SAMPLES_ERCC_MAPPING_DIR} 

#************************************* All samples, WT, TG TOFU

 # tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir>
tofu ${SAMPLES} ${SAMPLES_CLUSTER_DIR} ${SAMPLES_MM10_MAPPING_DIR} ${SAMPLES_MM10_TOFU_DIR}
tofu ${SAMPLES} ${SAMPLES_CLUSTER_DIR} ${SAMPLES_ERCC_MAPPING_DIR} ${SAMPLES_ERCC_TOFU_DIR}

#************************************* Merge RNASeq fastq files
# RNASeq_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name> <samples...>
TG_SAMPLES=(K24 L22 M20 O24 P22 Q20 S24	T22	K22	L20	M18	O22	P20	Q18	S22	T20	K20	L18	M24	O20	P18	Q24	S20	T18	K18	L24	M22	O18	P24	Q22	S18	T24)
WT_SAMPLES=(K17 L21 M19 K23 P21 Q19 M21 T21 K21 L19 M17 O21 P19 Q17 S21 T19 K19 L17 M23 O19 P17 Q23 S19 T17 O23 L23 Q21 O17 P23 S23 S17 T23)
ALL_SAMPLES+=( "${TG_SAMPLES[@]}" "${WT_SAMPLES[@]}" )

RNASeq_merge_fastq $RNASeq_Filtered $KALLISTO All_RNASeq ${ALL_SAMPLES[@]}
RNASeq_merge_fastq $RNASeq_Filtered $KALLISTO TG_RNASeq ${TG_SAMPLES[@]}
RNASeq_merge_fastq $RNASeq_Filtered $KALLISTO WT_RNASeq ${WT_SAMPLES[@]}


run_kallisto(){
    source activate sqanti2

    echo "Processing Kallisto for $1"
    cd $4
    kallisto version
    time kallisto index -i $1_Kallisto.idx $2 2> $1_Kallisto.index.log
    time kallisto quant -i $4/$1_Kallisto.idx --fr-stranded $3/$1_R1.fq --rf-stranded $3/$1_R2.fq -o $4 2> $1_Kallisto.quant.log
    mv abundance.tsv $1".abundance.tsv"

    # problem: retained co-ordinates, which does not input well into SQANTI2
    echo "Kallisto original $1.abundance.tsv"
    head $1".abundance.tsv"
    # solution: retain the PB.ID
    while read line ; do
      first=$( echo "$line" |cut -d\| -f1 ) # each read line, and remove the first part i.e. PacBio ID
      rest=$( echo "$line" | cut -d$'\t' -f2-5 ) #save the remaining columns
      echo $first $rest # concatenate
    done < $4/$1".abundance.tsv" > $4/$1.temp.abundance.tsv

    header=$(head -n 1 $4/$1.abundance.tsv)
    sed -i '1d' $4/$1.temp.abundance.tsv # remove header of temp.file to be replaced
    echo $header > foo
    cat foo $4/$1.temp.abundance.tsv > $4/$1.mod.abundance.tsv

    echo "Kallisto $1.mod.abundance.tsv"
    head $4/$1.mod.abundance.tsv
    rm $1.temp.abundance.tsv
    rm foo

    source deactivate
}

# run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto All_RNASeq $All_Isoseq3_WKD/TOFU/mm10/All_Merged.collapsed.filtered.rep.fa $KALLISTO $KALLISTO
run_kallisto TG_RNASeq $TG_Isoseq3_WKD/TOFU/mm10/TG_Merged.collapsed.filtered.rep.fa $KALLISTO $KALLISTO
run_kallisto WT_RNASeq $WT_Isoseq3_WKD/TOFU/mm10/WT_Merged.collapsed.filtered.rep.fa $KALLISTO $KALLISTO


kallisto quant -i $KALLISTO/TG_RNASeq_Kallisto.idx --fr-stranded TG_RNASeq_R1.fq --rf-stranded TG_RNASeq_R2.fq -o $KALLISTO 2> $1_Kallisto.quant.log

#************************************* All samples, WT, TG Chain and Individual samples
# Prepare files for chaining as according to Liz's chain.py
cd $CHAIN; mkdir ERCC mm10; cd mm10; mkdir All_Merged TG_Merged WT_Merged
cd $CHAIN/ERCC; mkdir All_Merged TG_Merged WT_Merged

SAMPLE_NAMES=(All_Merged TG_Merged WT_Merged)
SAMPLES_MM10_CHAIN_INPUT_DIR=($CHAIN/mm10/All_Merged $CHAIN/mm10/TG_Merged $CHAIN/mm10/WT_Merged)
SAMPLES_ERCC_CHAIN_INPUT_DIR=($CHAIN/ERCC/All_Merged $CHAIN/ERCC/TG_Merged $CHAIN/ERCC/WT_Merged)
SAMPLES_MM10_TOFU_DIR=($All_Isoseq3_WKD/TOFU/mm10 $TG_Isoseq3_WKD/TOFU/mm10 $WT_Isoseq3_WKD/TOFU/mm10)
SAMPLES_ERCC_TOFU_DIR=($All_Isoseq3_WKD/TOFU/ERCC $TG_Isoseq3_WKD/TOFU/ERCC $WT_Isoseq3_WKD/TOFU/ERCC)

# make sure not in conda environment when renaming
for i in {0..2}; do
	sample=${SAMPLE_NAMES[$i]}
	echo "Processing $sample for mm10 Chaining"
	cd ${SAMPLES_MM10_CHAIN_INPUT_DIR[$i]}
	cp ${SAMPLES_MM10_TOFU_DIR[$i]}/{*$sample.collapsed.group.txt*,*$sample.collapsed.filtered.gff*,*$sample.collapsed.filtered.abundance.txt*,*$sample.collapsed.filtered.rep.fq*} .
	echo "Copied over files"
	ls
	rename $sample Sample *

	echo "Processing $sample for ERCC Chaining"
	cd ${SAMPLES_ERCC_CHAIN_INPUT_DIR[$i]}
	cp ${SAMPLES_ERCC_TOFU_DIR[$i]}/{*$sample.collapsed.group.txt*,*$sample.collapsed.filtered.gff*,*$sample.collapsed.filtered.abundance.txt*,*$sample.collapsed.filtered.rep.fq*} .
	echo "Copied over files"
	ls
	rename $sample Sample *
done


# Prepare configuration file
cat << EOF > $CHAIN/mm10/mm10_Chained_Configuration.txt
SAMPLE=All_Merged;$CHAIN/mm10/All_Merged
SAMPLE=TG_Merged;$CHAIN/mm10/TG_Merged
SAMPLE=WT_Merged;$CHAIN/mm10/WT_Merged

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF

cat << EOF > $CHAIN/ERCC/ERCC_Chained_Configuration.txt
SAMPLE=All_Merged;$CHAIN/mm10/All_Merged
SAMPLE=TG_Merged;$CHAIN/mm10/TG_Merged
SAMPLE=WT_Merged;$CHAIN/mm10/WT_Merged

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF

module load Miniconda2
source activate sqanti2_py3
echo "Processing samples for chain_samples.py"
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/cupcake/tofu/counting
cd $CHAIN/mm10/;python $CUPCAKE/chain_samples_py3.py $CHAIN/mm10/mm10_Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> mm10_Chained_Configuration.log
cd $CHAIN/ERCC/;python chain_samples.py $CHAIN/ERCC/ERCC_Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> ERCC_Chained_Configuration.log