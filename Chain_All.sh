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

# 23/10/2020: CHAIN ALL_MERGED WT_MERGED TG_MERGED, KALLISTO

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

################## chain individual samples 
# Prepare files for chaining as according to Liz's chain.py 

SAMPLE_NAMES=(O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23)
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/TOFU
cd $CHAIN_INDIVIDUAL; mkdir ${SAMPLE_NAMES[@]}

# make sure not in conda environment when renaming
for sample in ${SAMPLE_NAMES[@]}; do	
	echo "Processing $sample for mm10 Chaining"
	cd $CHAIN_INDIVIDUAL/$sample
	cp $TOFU/{*$sample.collapsed.group.txt*,*$sample.collapsed.filtered.gff*,*$sample.collapsed.filtered.abundance.txt*,*$sample.collapsed.filtered.rep.fq*} .
	echo "Copied over files"
	ls
	rename $sample Sample *	
done


# Prepare configuration file 
cat << EOF > $CHAIN_INDIVIDUAL/Chained_Configuration.txt 
SAMPLE=O18;$CHAIN_INDIVIDUAL/O18
SAMPLE=K18;$CHAIN_INDIVIDUAL/K18
SAMPLE=S18;$CHAIN_INDIVIDUAL/S18
SAMPLE=L22;$CHAIN_INDIVIDUAL/L22
SAMPLE=Q20;$CHAIN_INDIVIDUAL/Q20
SAMPLE=K24;$CHAIN_INDIVIDUAL/K24
SAMPLE=Q21;$CHAIN_INDIVIDUAL/Q21
SAMPLE=K17;$CHAIN_INDIVIDUAL/K17
SAMPLE=M21;$CHAIN_INDIVIDUAL/M21
SAMPLE=O23;$CHAIN_INDIVIDUAL/O23
SAMPLE=S23;$CHAIN_INDIVIDUAL/S23
SAMPLE=K23;$CHAIN_INDIVIDUAL/K23

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF

echo "Processing samples for chain_samples.py"
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/cupcake/tofu/counting
cd $CHAIN_INDIVIDUAL;python $CUPCAKE/chain_samples_py3.py $CHAIN_INDIVIDUAL/Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> Chained_Configuration.log

SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
python $SQANTI2_DIR/sqanti_qc2.py --version
python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf ../all_samples.chained.gff $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --expression ../all_samples.chained_count.txt

SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence



#************************************* RAREFACTION
DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/
TOFU_INPUT_DIR=($DIR/All_Merged/TOFU/mm10 $DIR/TG_Merged/TOFU/mm10 $DIR/WT_Merged/TOFU/mm10)
RAREFACTION_OUTPUT_DIR=($DIR/All_Merged/RAREFACTION $DIR/TG_Merged/RAREFACTION $DIR/WT_Merged/RAREFACTION)
SQANTI_INPUT_DIR=($DIR/All_Merged/SQANTI $DIR/TG_Merged/SQANTI $DIR/WT_Merged/SQANTI)
SAMPLES=(All_Merged TG_Merged WT_Merged)

# make_file_for_rarefaction <sample_name_prefix> <input_tofu_directory> <working_directory> <input_sqanti_directory> 
make_file_for_rarefaction(){
    CUPCAKE_ANNOTATION=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/annotation
	cd $3
	echo "Working with $1" 
	# make_file_for_subsampling_from_collapsed.py <sample_name_prefix>.input.file <sample_name_prefix>.output.file <sample_name_prefix>.classification.txt
	python $CUPCAKE_ANNOTATION/make_file_for_subsampling_from_collapsed.py -i $2/$1.collapsed.filtered -o $1.subsampling -m2 $4/$1.collapsed.filtered.rep_classification.txt
	python $CUPCAKE_ANNOTATION/subsample_with_category.py --by refisoform --min_fl_count 2 --step 1000 $1.subsampling.all.txt > $1.rarefaction.by_refisoform.min_fl_2.by_category.txt
}

make_file_for_rarefaction $SAMPLES $TOFU_INPUT_DIR $RAREFACTION_OUTPUT_DIR $SQANTI_INPUT_DIR

