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
#SBATCH --output=cupcake85.o
#SBATCH --error=cupcake85.e

# 15/12/2020: Merge all Samples and then demultiplex later for count matrix
# 06/01/2020: Rerun with reduced coverage for cupcake parameters and TAMA filte

#************************************* DEFINE GLOBAL VARIABLES
All_Merged=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Merged
CCS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/CCS
SAMPLES_NAMES=(Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)
RAW_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Rawdata
SAMPLES_LIST=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome
FASTA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/primer.fasta

RAW_DATA_SUBREADS=(
#1.Q21_Tg4510_WT_2mos
$RAW_DIR/m54082_180607_173058.subreads.bam
#2.O18_Tg4510_TG_2mos
$RAW_DIR/m54082_180605_141944.subreads.bam
#3.L22_Tg4510_TG_8months
$RAW_DIR/m54082_190306_083150.subreads.bam
#4.K18_Tg4510_TG_2months
$RAW_DIR/m54082_190307_045507.subreads.bam
#5.O23_Tg4510_WT_8mos
$RAW_DIR/m54082_190401_165425.subreads.bam
#6.S23_Tg4510_WT_8mos
$RAW_DIR/m54082_190403_135102.subreads.bam
#7.S18_Tg4510_TG_2mos
$RAW_DIR/m54082_190404_101400.subreads.bam
#8.K17_Tg4510_WT_2mos
$RAW_DIR/m54082_190405_063832.subreads.bam
#9.M21_Tg4510_WT_2mos
$RAW_DIR/m54082_190430_163756.subreads.bam
#10.K23_Tg4510_WT_8mos
$RAW_DIR/m54082_190524_145911.subreads.bam
#11.Q20_Tg4510_TG_8mos
$RAW_DIR/m54082_190527_173356.subreads.bam
#12.K24_Tg4510_TG_8mos
$RAW_DIR/m54082_190529_082942.subreads.bam
)

module load Miniconda2/4.3.21
#************************************* All samples, WT, TG merged at refine
# run_isoseq3_merge <sample_prefix_output_name> <working_directory> <input_ccs_file>
# Aim: Pipeline for Isoseq3.1.2 (Lima, Refine, Polish) on merged samples
# Prerequisite:
    # DEFINE SAMPLES_NAMES list in working script
    # DEFINE RAW_DATA_SUBREADS in working script (path directories of subreadset.xml files)
    # DEFINE FASTA in working script (primer fasta)
run_isoseq3_merge(){
    source activate isoseq3
    ccs --version
    lima --version
    isoseq3 --version

    echo "Processing isoseq3 merge of Sample $1"

    cd $2
    # Define variable "Merge_Samples" as a list of all samples, in order to find the specified ccs.bam (for dataset create ConsensusReadSet)
    Merge_Samples=$(echo ${SAMPLES_NAMES[@]})

    # Define variable "all_ccs_bams" for merged path-directory of all ccs samples (for dataset create ConsensusReadSet)
    echo "Merging ccs of samples $Merge_Samples"
    all_ccs_bams=$(
        for i in ${Merge_Samples[@]}; do
            ccs_bam_name=$(find $3 -name "*.ccs.bam" -exec basename \{} \; | grep ^$i )
            ccs_bam=$(find $3 -name "*.ccs.bam" | grep "$ccs_bam_name" )
            echo $ccs_bam
        done
    )
    echo "Merging the following ccs.bams"
    echo ${all_ccs_bams[@]}
    # dataset create --type ConsensusReadSet <output.ccs.consensusreadset.xml> <input_ccs_bams>
    dataset create --type ConsensusReadSet $1_ccs.consensusreadset.xml $all_ccs_bams

    echo "Processing LIMA for sample $1"
    # lima <input.ccs.merged.consensusreadset.xml> <input.primerfasta> <output.fl.bam>
    lima $1_ccs.consensusreadset.xml $FASTA $1.fl.bam --isoseq --dump-clips

    echo "Processing Isoseq3 refine, cluster and polish"
    # refine --require-polya <input..fl.primer_5p--primer_3p.bam> <input.primer.fasta> <output.flnc.bam>
    isoseq3 refine --require-polya $1.fl.primer_5p--primer_3p.bam $FASTA $1.flnc.bam

    # cluster <input.flnc.bam> <output.unpolished.bam>
    isoseq3 cluster $1.flnc.bam $1.clustered.bam -j 32 --verbose --use-qvs 2> $1.cluster.log
    gunzip *.gz

    source deactivate
}

# run_post_isoseq3 <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
# Prerequisite: mm10 cage peak
run_post_isoseq3(){

    CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake
    SQANTI2_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    ANNOTATION=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/annotation
    SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence

    source activate sqanti2_py3

    # convert between fasta and fastq for downstream process
    echo "fasta2fastq conversion"
    #python $SEQUENCE/fa2fq.py $2/$1.clustered.hq.fasta

    samtools --version # echo version
    echo "Minimap2 version:" $(minimap2 --version) # echo version

    echo "Processing Sample $1 for Minimap2 and sort"
    #cd $3 #cd to $MAP directory for output
    #minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $REFERENCE/mm10.fa $2/$1.clustered.hq.fastq > $1.sam 2> $1.map.log
    #samtools sort -O SAM $1.sam > $1.sorted.sam

    head $CUPCAKE/README.md
    echo "Processing Sample $1 for ToFU collapse"
    cd $4 #cd to $TOFU directory for output
    if python $CUPCAKE/cupcake/tofu/collapse_isoforms_by_sam.py --input $2/$1.clustered.hq.fastq --fq -s $3/$1.sorted.sam --dun-merge-5-shorter -o $1 2> $1.collapse.log; then
        echo "Processing Sample $1 for collapse successful"
    else
        echo "Processing Sample $1 for collapse failed"
    fi

    if python $CUPCAKE/cupcake/tofu/get_abundance_post_collapse.py $1.collapsed $2/$1.clustered.cluster_report.csv 2> $1.abundance.log; then
        echo "Processing Sample $1 for getting abundance successful"
    else
        echo "Processing Sample $1 for getting abundance failed"
    fi

    if python $CUPCAKE/cupcake/tofu/filter_away_subset.py $1.collapsed 2> $1.filter.log; then
        echo "Processing Sample $1 for filtering successful"
    else
        echo "Processing Sample $1 for filtering failed"
    fi

    seqtk seq -a $1.collapsed.filtered.rep.fq > $1.collapsed.filtered.rep.fa
    source deactivate
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
    mv *png* Unfiltered_PNG/

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

# tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir> <default/reduced>
tofu_params(){
    source activate cupcake
    cd $4

    # Collapse
    if [ $5 = "default" ]; then
      echo "Processing Sample $1 for TOFU, with coverage 99% and identity 95%"
      collapse_isoforms_by_sam.py -c 0.99 -i 0.95 --input $2/$1.clustered.hq.fastq --fq -s $3/$1.sorted.sam --dun-merge-5-shorter -o $1 &> $1.collapse.log
    elif [ $5 = "reduced" ]; then
	    echo "Processing Sample $1 for TOFU, with coverage 95% and identity 95%"
      collapse_isoforms_by_sam.py -c 0.95 -i 0.95 --input $2/$1.clustered.hq.fastq --fq -s $3/$1.sorted.sam --dun-merge-5-shorter -o $1 &> $1.collapse.log
    elif [ $5 = "sigreduced" ]; then
	    echo "Processing Sample $1 for TOFU, with coverage 85% and identity 95%"
      collapse_isoforms_by_sam.py -c 0.85 -i 0.95 --input $2/$1.clustered.hq.fastq --fq -s $3/$1.sorted.sam --dun-merge-5-shorter -o $1 &> $1.collapse.log
	  else
        echo "5th parameter has to be set to default or reduced"
    fi

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


# demux <input path read.stat file> <input path of samples file> <path of output>
# run Cupcake_Demultiplex.R, read in read.stat file from cupcake collapse output and count abundance of each sample based on CCS_ID
demux(){
    source activate sqanti2_py3
    DEMUXFUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

    # Rscript.R <input path read.stat file> <input path of samples file> <path of output>
    Rscript $DEMUXFUNCTIONS/Cupcake_Demultiplex.R $1 $2 $3
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

# TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_name> <output_dir>
TAMA_sqanti_filter(){
  source activate sqanti2_py3
  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
  # Rscript .R <path/tama_filtered_output> <sqanti_filtered_dir> <sqanti_output_name> <output_dir>
  Rscript $GENERALFUNC/tama_sqanti_classgtfsubset.R $1 $2 $3 $4
  # extract fasta sequence based on the pbid (https://www.biostars.org/p/319099/)
  # script.py <path/sqanti_filtered.fasta> <path/retained_pbid_tama.txt> <path/output.fasta>
  cd $4
  awk '{ print $4 }' $1| cut -d ";" -f 2  > tama_retained_pbid.txt
  python $GENERALFUNC/tama_sqanti_fastasubset.py $2/$3".filtered_classification.filtered_lite.fasta" $4/tama_retained_pbid.txt $4/$3"_sqantifiltered_tamafiltered_classification.fasta"
}

#************************************* Rerun with changed cupcake parameters to 85% (08/01/2021)
cd $All_Merged; mkdir DEMUX_c85; cd DEMUX_c85; mkdir TOFU TAMA_Filter SQANTI
# tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir> <default/reduced/sigreduced>
# demux <input path read.stat file> <input path of samples file> <path of output>
# run_sqanti2 <sample_name> <working_directory> <collapsed.gtf> <count.txt>
# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
#tofu_params WholeIsoSeq $All_Merged/IsoSeq $All_Merged/MAPPING $All_Merged/DEMUX_c85/TOFU sigreduced
demux $All_Merged/DEMUX_c85/TOFU/WholeIsoSeq.collapsed.read_stat.txt $SAMPLES_LIST/Demultiplex_IsoSeq_Whole.csv $All_Merged/DEMUX_c85/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
run_sqanti2 WholeIsoSeq.collapsed.filtered $All_Merged/DEMUX_c85/SQANTI $All_Merged/DEMUX_c85/TOFU/WholeIsoSeq.collapsed.filtered.gff $All_Merged/DEMUX_c85/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
TAMA_remove_fragments $All_Merged/DEMUX_c85/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite.gtf WholeIsoSeq $All_Merged/DEMUX_c85/TAMA_Filter