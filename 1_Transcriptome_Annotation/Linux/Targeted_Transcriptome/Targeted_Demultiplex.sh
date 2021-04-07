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
#SBATCH --output=Targeted_Demultiplex.o
#SBATCH --error=Targeted_Demultiplex.e

# 15/12/2020: Merge all Samples and then demultiplex later for count matrix

############# Alternatively cat multiple refine files and therefore no need to cluster (which takes a lot of time)
#************************************* DEFINE GLOBAL VARIABLES
# File directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
DEMUX_SCRIPT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome/
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/cupcake/tofu/counting
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/IsoSeq/REFINE
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged/CLUSTER
All_Samples=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/All_Targeted_Merged
Targeted_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/
GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/
cd $All_Samples; mkdir CLUSTER MAPPING TOFU SQANTI

#************************************* All samples, WT, TG merged at refine
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
source $FUNCTIONS/Post_Isoseq3_Function.sh
module load Miniconda2
module load R

#run_pipeline_to_chain <Sample_Name> <List.of.Samples...>
run_pipeline_to_tofu(){

  SAMPLES=$(echo "${@:2}")
  echo "Processing together: $SAMPLES"

  # merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
  # convert_fa2fq <file_name> <input_dir>
  # run_minimap2 <prefix_sample> <input_dir> <ERCC/mm10> <output_dir>
  # tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir>
  # chain_prepare <ERCC/MOUSE>
  #merging_at_refine $REFINE $All_Samples/CLUSTER $1 ${SAMPLES[@]}
  #convert_fa2fq $1".clustered.hq.fasta" $All_Samples/CLUSTER
  #run_minimap2 $1 $All_Samples/CLUSTER mm10 $All_Samples/MAPPING
  #on_target_rate $Targeted_dir/Probes/FINAL_MOUSE.bed $All_Samples/CLUSTER/All_Targeted_Merged.clustered.hq.fasta $All_Samples/MAPPING/All_Targeted_Merged.clustered.hq.fastq.sam $All_Samples/MAPPING/All_Targeted_Merged.fasta.sam.probe_hit.txt
  #tofu $1 $All_Samples/CLUSTER $All_Samples/MAPPING $All_Samples/TOFU

  # demultiplex collapsed.read_stat.txt
  Rscript $DEMUX_SCRIPT/Demultiplex_Cupcake.R
}

# run_sqanti2 <sample_name> <working_directory> <collapsed.gtf> <count.txt>
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

    mkdir Unfiltered_PNG
    mv *png* Unfiltered_PNG

    if python $SQANTI2_DIR/sqanti_filter2.py \
        $1"_classification.txt" \
        $1"_corrected.fasta"  \
        $1"_corrected.gtf" \
        -a 0.6 -c 3 2> $1.sqanti.filter.log; then
        echo "Processing Sample $1 for SQANTI2 filter successful"
    else
        echo "Processing Sample $1 for SQANTI2 filter failed"
    fi

    mkdir filtered_PNG
    mv *png* filtered_PNG/


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
    seqtk seq -a $1".collapsed.rep.fq" > $1".collapsed.rep.fa"
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


run_pipeline_to_tofu All_Targeted_Merged K17 K18 K19 K20 K21 K23 K24 L18 L22 M21 O18 O22 O23 P19 Q17 Q18 Q20 Q21 Q23 S18 S19 S23 T18 T20

# run_sqanti2_QC <sample_name> <working_directory> <collapsed.gtf> <count.txt>
# TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_dir>
# Rscript Targeted_Cupcake_Demultiplex.R <path/refine_dir> <cluster_report> <path/collapsed_read_stat file> <output_file>
#Rscript $GENERALFUNC/Targeted_Cupcake_Demultiplex.R $REFINE $CLUSTER/All_Targeted_Merged.clustered.cluster_report.csv $All_Samples/TOFU/All_Targeted_Merged.collapsed.read_stat.txt  $All_Samples/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt
run_sqanti2 All_Targeted_Merged.collapsed.filtered $All_Samples/SQANTI $All_Samples/TOFU/All_Targeted_Merged.collapsed.filtered.gff $All_Samples/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt
# run Rscript: Demultiplex_Cupcake_Subset.R
TAMA_remove_fragments $All_Samples/SQANTI/All_Targeted_Merged.collapsed.filtered_classification.filtered_lite.gtf All_Targeted_Merged_Subset $All_Samples/TAMA
TAMA_sqanti_filter $All_Samples/TAMA/All_Targeted_Merged.bed $All_Samples/SQANTI All_Targeted_Merged.collapsed.filtered_classification.filtered_lite_classification.txt All_Targeted_Merged.collapsed.filtered_classification.filtered_lite.gtf All_Targeted_Merged.collapsed.filtered_classification.filtered_lite.fasta All_Targeted_Merged $All_Samples/TAMA_SQANTI_FILTER

### run from mapping; alternative pipeline without cupcake filtering
#mkdir $All_Samples/Alternative_Pipeline; cd $All_Samples/Alternative_Pipeline; mkdir TOFU SQANTI TAMA SQANTI_TAMA_FILTER
tofu_collapse All_Targeted_Merged $All_Samples/CLUSTER $All_Samples/MAPPING $All_Samples/Alternative_Pipeline/TOFU 0.99 0.95
Rscript $GENERALFUNC/Targeted_Cupcake_Demultiplex.R $REFINE $CLUSTER/All_Targeted_Merged.clustered.cluster_report.csv $All_Samples/Alternative_Pipeline/TOFU/All_Targeted_Merged.collapsed.read_stat.txt  $All_Samples/Alternative_Pipeline/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt
run_sqanti2 All_Targeted_Merged.collapsed $All_Samples/Alternative_Pipeline/SQANTI $All_Samples/Alternative_Pipeline/TOFU/All_Targeted_Merged.collapsed.gff $All_Samples/Alternative_Pipeline/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt
TAMA_remove_fragments $All_Samples/Alternative_Pipeline/SQANTI/All_Targeted_Merged.collapsed_classification.filtered_lite.gtf All_Targeted_Merged $All_Samples/Alternative_Pipeline/TAMA
TAMA_sqanti_filter $All_Samples/Alternative_Pipeline/TAMA/All_Targeted_Merged.bed $All_Samples/Alternative_Pipeline/SQANTI All_Targeted_Merged.collapsed_classification.filtered_lite_classification.txt All_Targeted_Merged.collapsed_classification.filtered_lite.gtf All_Targeted_Merged.collapsed_classification.filtered_lite.fasta All_Targeted_Merged $All_Samples/Alternative_Pipeline/SQANTI_TAMA_FILTER

### run from mapping; alternative pipeline with collapse threshold at 85% and cupcake filtering
mkdir $All_Samples/Alternative_Pipeline/c85; cd $All_Samples/Alternative_Pipeline/c85; mkdir TOFU SQANTI TAMA SQANTI_TAMA_FILTER
tofu_collapse All_Targeted_Merged $All_Samples/CLUSTER $All_Samples/MAPPING $All_Samples/Alternative_Pipeline/c85/TOFU 0.85 0.95
Rscript $GENERALFUNC/Targeted_Cupcake_Demultiplex.R $REFINE $CLUSTER/All_Targeted_Merged.clustered.cluster_report.csv $All_Samples/Alternative_Pipeline/c85/TOFU/All_Targeted_Merged.collapsed.read_stat.txt  $All_Samples/Alternative_Pipeline/c85/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt
run_sqanti2 All_Targeted_Merged.collapsed $All_Samples/Alternative_Pipeline/c85/SQANTI $All_Samples/Alternative_Pipeline/c85/TOFU/All_Targeted_Merged.collapsed.gff $All_Samples/Alternative_Pipeline/c85/TOFU/All_Targeted_Merged.Demultipled_Abundance.txt
TAMA_remove_fragments $All_Samples/Alternative_Pipeline/c85/SQANTI/All_Targeted_Merged.collapsed_classification.filtered_lite.gtf All_Targeted_Merged $All_Samples/Alternative_Pipeline/c85/TAMA
TAMA_sqanti_filter $All_Samples/Alternative_Pipeline/c85/TAMA/All_Targeted_Merged.bed $All_Samples/Alternative_Pipeline/c85/SQANTI All_Targeted_Merged.collapsed_classification.filtered_lite_classification.txt All_Targeted_Merged.collapsed_classification.filtered_lite.gtf All_Targeted_Merged.collapsed_classification.filtered_lite.fasta All_Targeted_Merged $All_Samples/Alternative_Pipeline/c85/SQANTI_TAMA_FILTER
