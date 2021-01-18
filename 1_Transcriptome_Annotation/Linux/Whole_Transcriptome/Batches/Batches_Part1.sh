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
run_sqanti2 All_Merged $SQANTI $TOFU/All_Merged.collapsed.filtered.gff $TOFU/All_Merged_demultiplexed_abundance.csv

# ERCC
SQANTI=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/SQANTIDEMUX
SQANTIAll=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/SQANTIAll
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/TOFU
run_sqanti2 All_Merged.collapsed.filtered $SQANTI $TOFU/All_Merged.collapsed.filtered.gff $SQANTI/Mouse.Demultiplexed_Abundance.txt ERCC
run_sqanti2 All_Merged.collapsed.filtered $SQANTIAll $TOFU/All_Merged.collapsed.filtered.gff $TOFU/All_Merged.collapsed.filtered.abundance.txt ERCC

#2b. Tama filter to remove redundant transcripts
# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
cd $All_Samples/ERCC/TOFU;mkdir TAMA_Filter
TAMA_remove_fragments $All_Samples/ERCC/TOFU/All_Merged.collapsed.filtered.gff All_Merged $All_Samples/ERCC/TOFU/TAMA_Filter
TAMA_remove_fragments $All_Samples/ERCC/SQANTIDEMUX/All_Merged.collapsed.filtered_classification.filtered_lite.gtf All_Merged_postsqanti $All_Samples/ERCC/TOFU/TAMA_Filter

# some ERCCS in minimap2 output but not CUPCAKE
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/CLUSTER
REFERENCE_ERCC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/ERCC
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/MAPPING
minimap2 -t 32 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $REFERENCE_ERCC/ERCC92.fa $CLUSTER/All_Merged.clustered.hq.fastq > All_Merged.clustered.hq.fastq.paf 2> All_Merged.clustered.hq.fastq.paf.log

htsbox samview -pS All_Merged.clustered.hq.fastq.sam > All_Merged.clustered.hq.fastq.paf
awk -F'\t' '{if ($6=="*") {print $0}}' All_Merged.clustered.hq.fastq.paf > All_Merged.notread.clustered.hq.fastq.paf
awk -F'\t' '{if ($6!="*") {print $0}}' All_Merged.clustered.hq.fastq.paf > All_Merged.filtered.clustered.hq.fastq.paf
awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' All_Merged.filtered.clustered.hq.fastq.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > All_Merged_reads_with_alignment_statistics.txt

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

cd $All_Samples/ERCC;mkdir TOFU_ADJUST
cd TOFU_ADJUST; mkdir SQANTI
DEMUX_ERCC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/SQANTIDEMUX
tofu_params All_Merged $All_Samples/CLUSTER $All_Samples/ERCC/MAPPING $All_Samples/ERCC/TOFU_ADJUST
run_sqanti2 All_Merged $All_Samples/ERCC/TOFU_ADJUST/SQANTI $All_Samples/ERCC/TOFU_ADJUST/All_Merged.collapsed.filtered.gff $DEMUX_ERCC/All_Merged_demultiplexed_abundance.csv

cd $All_Samples/ERCC/TOFU_ADJUST;mkdir TAMA_Filter
TAMA_remove_fragments $All_Samples/ERCC/TOFU_ADJUST/SQANTI/All_Merged.collapsed.filtered_classification.filtered_lite.gtf All_Merged_postsqanti $All_Samples/ERCC/TOFU_ADJUST/TAMA_Filter