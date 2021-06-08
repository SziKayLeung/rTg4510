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
#SBATCH --output=Batches_Part1.o
#SBATCH --error=Batches_Part1.e

# 14/12/2020: Merging All, TG and WT samples and processing to CHAIN
# 04/01/2021: Rerunning chaining with latest v18.1 cupcake
# 22/01/2021: Rerun ERCC with new pipeline (skip cupcake filter, include TAMA filter)

#************************************* DEFINE GLOBAL VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/cupcake/tofu/counting
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/REFINE
All_Samples=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples
#cd $All_Samples; mkdir CLUSTER ERCC MOUSE
#cd $All_Samples/ERCC; mkdir MAPPING TOFU CHAIN SQANTI
#cd $All_Samples/MOUSE; mkdir MAPPING TOFU CHAIN SQANTI

#************************************* All samples, WT, TG merged at refine
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
source $FUNCTIONS/Post_Isoseq3_Function.sh

module load Miniconda2
# chain_prepare <ERCC/MOUSE>
chain_prepare_ERCC(){
  source deactivate

  ALL_SAMPLES_NAMES=(All_Merged TG_Merged WT_Merged)
  cd $All_Samples/ERCC/CHAIN; mkdir All_Merged TG_Merged WT_Merged

  for sample in ${ALL_SAMPLES_NAMES[@]}; do
    echo "Processing $sample for ERCC Chaining"
  	cd $All_Samples/ERCC/CHAIN/$sample; cp $All_Samples/ERCC/TOFU/{$sample".collapsed.group.txt",$sample".collapsed.filtered.gff",$sample".collapsed.filtered.abundance.txt",$sample".collapsed.filtered.rep.fq"} .
  	echo "Copied over files"; ls
  	rename $sample Sample *
  done

cat << EOF > $All_Samples/ERCC/CHAIN/Chained_Configuration.txt
SAMPLE=All_Merged;$All_Samples/ERCC/CHAIN/All_Merged
SAMPLE=TG_Merged;$All_Samples/ERCC/CHAIN/TG_Merged
SAMPLE=WT_Merged;$All_Samples/ERCC/CHAIN/WT_Merged

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF

  source activate cupcake

}

chain_prepare_MOUSE(){
  source deactivate

  ALL_SAMPLES_NAMES=(All_Merged TG_Merged WT_Merged)
  cd $All_Samples/MOUSE/CHAIN; mkdir All_Merged TG_Merged WT_Merged

  for sample in ${ALL_SAMPLES_NAMES[@]}; do
    echo "Processing $sample for MOUSE Chaining"
  	cd $All_Samples/MOUSE/CHAIN/$sample; cp $All_Samples/MOUSE/TOFU/{$sample".collapsed.group.txt",$sample".collapsed.filtered.gff",$sample".collapsed.filtered.abundance.txt",$sample".collapsed.filtered.rep.fq"} .
  	echo "Copied over files"; ls
  	rename $sample Sample *
  done

cat << EOF > $All_Samples/MOUSE/CHAIN/Chained_Configuration.txt
SAMPLE=All_Merged;$All_Samples/MOUSE/CHAIN/All_Merged
SAMPLE=TG_Merged;$All_Samples/MOUSE/CHAIN/TG_Merged
SAMPLE=WT_Merged;$All_Samples/MOUSE/CHAIN/WT_Merged

GROUP_FILENAME=Sample.collapsed.group.txt
GFF_FILENAME=Sample.collapsed.filtered.gff
COUNT_FILENAME=Sample.collapsed.filtered.abundance.txt
FASTQ_FILENAME=Sample.collapsed.filtered.rep.fq
EOF

  source activate cupcake

}
#run_pipeline_to_chain <Sample_Name> <List.of.Samples...>
run_pipeline_to_chain(){

  SAMPLES=$(echo "${@:2}")
  echo "Processing together: $SAMPLES"

  #export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/cDNA_Cupcake/cupcake/tofu/counting

  # merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
  # convert_fa2fq <file_name> <input_dir>
  # run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
  # tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir>
  # chain_prepare <ERCC/MOUSE>
  #merging_at_refine $REFINE $All_Samples/CLUSTER $1 ${SAMPLES[@]}
  #convert_fa2fq $1".clustered.hq.fasta" $All_Samples/CLUSTER

  ## ERCC
  #run_minimap2 $1 $All_Samples/CLUSTER ERCC $All_Samples/ERCC/MAPPING
  #tofu $1 $All_Samples/CLUSTER $All_Samples/ERCC/MAPPING $All_Samples/ERCC/TOFU
  #chain_prepare_ERCC
  source activate cupcake
  cd $All_Samples/ERCC/CHAIN; chain_samples.py $All_Samples/ERCC/CHAIN/Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> Chained_Configuration.log
  source deactivate

  ## MOUSE
  #run_minimap2 $1 $All_Samples/CLUSTER mm10 $All_Samples/MOUSE/MAPPING
  #tofu $1 $All_Samples/CLUSTER $All_Samples/MOUSE/MAPPING $All_Samples/MOUSE/TOFU
  #chain_prepare_MOUSE
  source activate cupcake
  cd $All_Samples/MOUSE/CHAIN; chain_samples.py $All_Samples/MOUSE/CHAIN/Chained_Configuration.txt count_fl --dun-merge-5-shorter 2> Chained_Configuration.log
  source deactivate
}

# run_sqanti2 <sample_name> <working_directory> <collapsed.gtf> <count.txt> <mm10/ERCC>
run_sqanti2(){

    source activate sqanti2_py3
    SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    CAGE_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CAGE

    echo "Processing Sample $1 for SQANTI2 QC"
    cd $2

    export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
    python $SQANTI2_DIR/sqanti_qc2.py -v

    if [ $5 == "mm10" ]; then
      if python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
          --gtf $3 \
          $REFERENCE/gencode.vM22.annotation.gtf  \
          $REFERENCE/mm10.fa \
          --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed \
          --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt \
          --fl_count $4 \
          2> $1.sqanti.qc.log; then
          echo "Processing Sample $1 for SQANTI2 successful"
      else
          echo "Processing Sample $1 for SQANTI2 failed"
      fi
     elif [ $5 == "ERCC" ]; then
       if python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
          --gtf $3 \
          $REFERENCE/ERCC/ERCC92.gtf  \
          $REFERENCE/ERCC/ERCC92.fa \
          --fl_count $4 \
          2> $1.sqanti.qc.log; then
          echo "Processing Sample $1 for SQANTI2 successful"
      else
          echo "Processing Sample $1 for SQANTI2 failed"
      fi
    else
        echo "ERROR: Require 5th argument to be either: mm10 or ERCC"
        return
    fi

    mkdir Unfiltered_PNG; mv *png* Unfiltered_PNG/

    if python $SQANTI2_DIR/sqanti_filter2.py \
        $1"_classification.txt" \
        $1"_corrected.fasta"  \
        $1"_corrected.gtf" \
        -a 0.6 -c 3 2> $1.sqanti.filter.log; then
        echo "Processing Sample $1 for SQANTI2 filter successful"
    else
        echo "Processing Sample $1 for SQANTI2 filter failed"
    fi

    mkdir filtered_PNG; mv *png* filtered_PNG/
    source deactivate

}

# TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
TAMA_sqanti_filter(){
  source activate sqanti2_py3
  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
  #Rscript .R <path/tama_filtered_output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <output_prefix_name> <output_dir>
  Rscript $GENERALFUNC/TAMA/tama_sqanti_classgtfsubset.R $1 $2 $3 $4 $6 $7
  # extract fasta sequence based on the pbid (https://www.biostars.org/p/319099/)
  # script.py <path/sqanti_filtered.fasta> <path/retained_pbid_tama.txt> <path/output.fasta>
  cd $7
  awk '{ print $4 }' $1| cut -d ";" -f 2  > tama_retained_pbid.txt
  python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $2/$5 $7/tama_retained_pbid.txt $7/$6"_sqantifiltered_tamafiltered_classification.fasta"
}



#1. chained
#run_pipeline_to_chain All_Merged O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23
#run_pipeline_to_chain TG_Merged O18 K18 S18 L22 Q20 K24
#run_pipeline_to_chain WT_Merged Q21 K17 M21 O23 S23 K23
#run_sqanti2 all_samples.chained $All_Samples/ERCC/SQANTI $All_Samples/ERCC/CHAIN/all_samples.chained.gff $All_Samples/ERCC/CHAIN/all_samples.chained_count.txt ERCC
#run_sqanti2 all_samples.chained $All_Samples/MOUSE/SQANTI $All_Samples/MOUSE/CHAIN/all_samples.chained.gff $All_Samples/MOUSE/CHAIN/all_samples.chained_count.txt mm10

#2. without chaining: demultiplex
# run_sqanti2 <sample_name> <working_directory> <collapsed.gtf> <count.txt>
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/TOFU
SQANTI=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/SQANTI_POSTTOFU
#run_sqanti2 All_Merged $SQANTI $TOFU/All_Merged.collapsed.filtered.gff $TOFU/All_Merged_demultiplexed_abundance.csv

# ERCC
SQANTI=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/SQANTIDEMUX
SQANTIAll=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/SQANTIAll
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/TOFU
#run_sqanti2 All_Merged.collapsed.filtered $SQANTI $TOFU/All_Merged.collapsed.filtered.gff $SQANTI/Mouse.Demultiplexed_Abundance.txt ERCC
#run_sqanti2 All_Merged.collapsed.filtered $SQANTIAll $TOFU/All_Merged.collapsed.filtered.gff $TOFU/All_Merged.collapsed.filtered.abundance.txt ERCC

#2b. Tama filter to remove redundant transcripts
# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
#cd $All_Samples/ERCC/TOFU;mkdir TAMA_Filter
#TAMA_remove_fragments $All_Samples/ERCC/SQANTIDEMUX/All_Merged.collapsed.filtered_classification.filtered_lite.gtf All_Merged_postsqanti $All_Samples/ERCC/TOFU/TAMA_Filter

# 2c. Filter SQANTI classification file after TAMA filter
# TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
#TAMA_sqanti_filter $All_Samples/ERCC/TOFU/TAMA_Filter/All_Merged_postsqanti.bed $All_Samples/ERCC/SQANTIDEMUX All_Merged.collapsed.filtered_classification.filtered_lite_classification.txt All_Merged.collapsed.filtered_classification.filtered_lite_classification.gtf All_Merged.collapsed.filtered_classification.filtered_lite_classification.fasta All_Merged $All_Samples/ERCC/TOFU/TAMA_SQANTI_FILTER


# some ERCCS in minimap2 output but not CUPCAKE
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/CLUSTER
REFERENCE_ERCC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/ERCC
#cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/MAPPING
#minimap2 -t 32 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $REFERENCE_ERCC/ERCC92.fa $CLUSTER/All_Merged.clustered.hq.fastq > All_Merged.clustered.hq.fastq.paf 2> All_Merged.clustered.hq.fastq.paf.log

#htsbox samview -pS All_Merged.clustered.hq.fastq.sam > All_Merged.clustered.hq.fastq.paf
#awk -F'\t' '{if ($6=="*") {print $0}}' All_Merged.clustered.hq.fastq.paf > All_Merged.notread.clustered.hq.fastq.paf
#awk -F'\t' '{if ($6!="*") {print $0}}' All_Merged.clustered.hq.fastq.paf > All_Merged.filtered.clustered.hq.fastq.paf
#awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' All_Merged.filtered.clustered.hq.fastq.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > All_Merged_reads_with_alignment_statistics.txt

#2d. cupcake collapse coverage (0.95) and identity (0.95 - default) readjustment
tofu_params(){
    source activate cupcake
    cd $4
    echo "Processing Sample $1 for TOFU, with coverage 95% and identity 95%"
    # Collapse
    collapse_isoforms_by_sam.py -c 0.95 -i 0.95 --input $2/$1.clustered.hq.fastq --fq -s $3/$1.clustered.hq.fastq.sorted.sam --dun-merge-5-shorter -o $1 \
    &> $1.collapse.log

    # Create Abundance Script of full-length transcripts
    get_abundance_post_collapse.py $1.collapsed $2/$1.clustered.cluster_report.csv 2> $1.abundance.log

    # Remove degraded isoforms (default setting)
    filter_away_subset.py $1.collapsed 2> $1.filter.log
    source activate

    source activate sqanti2_py3
    # convert rep.fq to rep.fa for SQANTI2 input
    seqtk seq -a $1.collapsed.filtered.rep.fq > $1.collapsed.filtered.rep.fa
    echo "Processing Sample $1 for TOFU successful"
    source deactivate
}

# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# remove short fragments from post tofu
# Prerequisite: Require TAMA_prepare.R to change column 4 of bed12 to gene_name: transcript_name for correct file TAMA format
TAMA_remove_fragments(){

    TAMAFUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/TAMA
    TAMA_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama/tama_go/filter_transcript_models
    source activate sqanti2
    cd $3

    # convert gtf to bed12
	  gtfToGenePred $1 $2.genepred
	  genePredToBed $2.genepred > $2.bed12
    awk -F'\t' '{print $1,$2,$3,$4,"40",$6,$7,$8,"255,0,0",$10,$11,$12}' $2.bed12| sed s/" "/"\t"/g|sed s/",\t"/"\t"/g|sed s/",$"/""/g > Tama_$2.bed12

    # Rscript script.R <name>.bed12 <input_dir>
	  Rscript $TAMAFUNCTIONS/TAMA_Merge_Prepare.R Tama_$2 $3
	  python $TAMA_DIR/tama_remove_fragment_models.py -f Tama_$2_mod.bed12 -o $2

	  rm *Tama*

    echo "Number of isoforms filtered by TAMA:"
    wc -l $2"_discarded.txt"

	  source deactivate
}

#cd $All_Samples/ERCC;mkdir TOFU_ADJUST
#cd TOFU_ADJUST; mkdir SQANTI
#DEMUX_ERCC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/SQANTIDEMUX
#tofu_params All_Merged $All_Samples/CLUSTER $All_Samples/ERCC/MAPPING $All_Samples/ERCC/TOFU_ADJUST
#run_sqanti2 All_Merged.collapsed $All_Samples/ERCC/TOFU_ADJUST/SQANTI $All_Samples/ERCC/TOFU_ADJUST/All_Merged.collapsed.filtered.gff $DEMUX_ERCC/All_Merged_demultiplexed_abundance.csv

#cd $All_Samples/ERCC/TOFU_ADJUST;mkdir TAMA_Filter
#TAMA_remove_fragments $All_Samples/ERCC/TOFU_ADJUST/SQANTI/All_Merged.collapsed.filtered_classification.filtered_lite.gtf All_Merged_postsqanti $All_Samples/ERCC/TOFU_ADJUST/TAMA_Filter

############################# 21/01/2021 New proposed pipeline with ERCCs
SAMPLES_LIST=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome
ERCC_ALT=$All_Samples/ERCC/Alternative_Pipeline
# demux <input path read.stat file> <input path of samples file> <path of output>
# run Cupcake_Demultiplex.R, read in read.stat file from cupcake collapse output and count abundance of each sample based on CCS_ID
demux(){
    source activate sqanti2_py3
    DEMUXFUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

    # Rscript.R <input path read.stat file> <input path of samples file> <path of output>
    Rscript $DEMUXFUNCTIONS/Cupcake_Demultiplex.R $1 $2 $3
}

# tofu_collapse <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir> <coverage_threshold> <length_threshold>
# <prefix_sample> = sample name for clustered.hq.fastq
# <input_dir> = directory path containing clustered.hq.fastq
tofu_collapse(){
    source activate cupcake
    cd $4
    echo "Processing Sample $1 for TOFU with coverage threshold at $5 and length threshold at $6"
    # Collapse
    collapse_isoforms_by_sam.py -c $5 -i $6 --input $2/$1.clustered.hq.fastq --fq -s $3/$1.clustered.hq.fastq.sorted.sam --dun-merge-5-shorter -o $1 &> $1.collapse.log

    # Create Abundance Script of full-length transcripts
    get_abundance_post_collapse.py $1.collapsed $2/$1.clustered.cluster_report.csv 2> $1.abundance.log
    source deactivate

    source activate sqanti2_py3
    # convert rep.fq to rep.fa for SQANTI2 input
    seqtk seq -a $1.collapsed.rep.fq > $1.collapsed.rep.fa
    echo "Processing Sample $1 for TOFU successful"
    source deactivate
}



# tofu_collapse <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir> <coverage_threshold> <length_threshold>
tofu(){
    source activate cupcake
    cd $4
    echo "Processing Sample $1 for TOFU with coverage threshold at $5 and length threshold at $6"
    # Collapse
    collapse_isoforms_by_sam.py -c $5 -i $6 --input $2/$1.clustered.hq.fastq --fq -s $3/$1.clustered.hq.fastq.sorted.sam --dun-merge-5-shorter -o $1 &> $1.collapse.log

    # Create Abundance Script of full-length transcripts
    get_abundance_post_collapse.py $1.collapsed $2/$1.clustered.cluster_report.csv 2> $1.abundance.log
    source deactivate

    # Remove degraded isoforms (default setting)
    filter_away_subset.py $1.collapsed 2> $1.filter.log
    source activate

    source activate sqanti2_py3
    # convert rep.fq to rep.fa for SQANTI2 input
    seqtk seq -a $1.collapsed.filtered.rep.fq > $1.collapsed.filtered.rep.fa
    echo "Processing Sample $1 for TOFU successful"
    source deactivate
}

#run_revised_pipeline <Sample_Name> <cluster_dir> <mapping_dir> <tofu_dir> <coverage_threshold> <sqanti_dir> <tama_dir> <santi_tama_dir> <List.of.Samples...>
run_revised_pipeline_ERCC(){

  # NOTE: error in SQANTI report due to no junctions file (there are all monoexonic); however, can continue as error only appears in generating sqanti report; classification is already complete
  SAMPLES=$(echo "${@:9}")
  echo "Processing together: $SAMPLES"

  #export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/cDNA_Cupcake/cupcake/tofu/counting

  # merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
  # convert_fa2fq <file_name> <input_dir>
  # run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
  # tofu_collapse <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir> <coverage_threshold> <length_threshold>
  # demux <input path read.stat file> <input path of samples file> <path of output>
  # TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
  # TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>

  #merging_at_refine $REFINE $All_Samples/CLUSTER $1 ${SAMPLES[@]}
  #convert_fa2fq $1".clustered.hq.fasta" $All_Samples/CLUSTER

  ## ERCC
  #run_minimap2 $1 $All_Samples/CLUSTER ERCC $All_Samples/ERCC/MAPPING
  # 95% coverage threshold
  #tofu_collapse $1 $2 $3 $4 $5 0.95
  tofu $1 $2 $3 $4 $5 0.95
  demux $4/$1".collapsed.read_stat.txt" $SAMPLES_LIST/Demultiplex_IsoSeq_Whole.csv $4/$1".Demultiplexed_Abundance.txt"
  run_sqanti2 $1.collapsed $ERCC_ALT/SQANTI/c95 $ERCC_ALT/TOFU/c95/$1".collapsed.gff" $4/$1".Demultiplexed_Abundance.txt" ERCC
  TAMA_remove_fragments $6/$1".collapsed_classification.filtered_lite.gtf" $1 $7
  TAMA_sqanti_filter $7/$1".bed" $6 $1".collapsed_classification.filtered_lite.txt" $1".collapsed_classification.filtered_lite.gtf" $1".collapsed_classification.filtered_lite.fasta" $1 $9
  echo "Pipeline run successfully"
}

# set up directories
#mkdir $ERCC_ALT
#dir=($ERCC_ALT/TOFU $ERCC_ALT/SQANTI $ERCC_ALT/TAMA_Filter $ERCC_ALT/SQANTI_TAMA_FILTER)
#for d in ${dir[@]};do mkdir $d; cd $d; mkdir c95 c99 c85; done

#run_revised_pipeline <Sample_Name> <cluster_dir> <mapping_dir> <tofu_dir> <coverage_threshold> <sqanti_dir> <tama_dir> <santi_tama_dir> <List.of.Samples...>
#run_revised_pipeline_ERCC All_Merged $All_Samples/CLUSTER $All_Samples/ERCC/MAPPING $ERCC_ALT/TOFU/c95 0.95 $ERCC_ALT/SQANTI/c95 $ERCC_ALT/TAMA_Filter/c95 $ERCC_ALT/SQANTI_TAMA_FILTER/c95 O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23
run_revised_pipeline_ERCC All_Merged $All_Samples/CLUSTER $All_Samples/ERCC/MAPPING $ERCC_ALT/TOFU/c99 0.99 $ERCC_ALT/SQANTI/c99 $ERCC_ALT/TAMA_Filter/c99 $ERCC_ALT/SQANTI_TAMA_FILTER/c99 O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23
#run_revised_pipeline_ERCC All_Merged $All_Samples/CLUSTER $All_Samples/ERCC/MAPPING $ERCC_ALT/TOFU/c85 0.85 $ERCC_ALT/SQANTI/c85 $ERCC_ALT/TAMA_Filter/c85 $ERCC_ALT/SQANTI_TAMA_FILTER/c85 O18 K18 S18 L22 Q20 K24 Q21 K17 M21 O23 S23 K23
