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
#SBATCH --output=All_Demultiplex.o
#SBATCH --error=All_Demultiplex.e

# 15/12/2020: Merge all Samples and then demultiplex later for count matrix
# 06/01/2020: Rerun with reduced coverage for cupcake parameters and TAMA filter
# 18/01/2021: Rerun with further reduced coverage to 85%

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
    if python $CUPCAKE/cupcake/tofu/collapse_isoforms_by_sam.py --input $2/$1.clustered.hq.fastq --fq -s $3/$1.sorted.sam --dun-merge-5-shorter -o $1 \
    2> $1.collapse.log; then
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
    Rscript $DEMUXFUNCTIONS/Cupcake_Demultiplex.R  $1 $2 $3
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

#************************************* Initial run (15/12/2021)
initial_run_cmd(){
run_isoseq3_merge WholeIsoSeq $All_Merged/IsoSeq $CCS
run_post_isoseq3 WholeIsoSeq $All_Merged/IsoSeq $All_Merged/MAPPING $All_Merged/TOFU
Demultiplex_Samples=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome/Demultiplex_IsoSeq_Whole.txt

cd $All_Merged/DEMUX
Rscript /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Whole_Transcriptome/ModifyClassifyReport.R $All_Merged/IsoSeq/WholeIsoSeq.flnc.report.csv $All_Merged/DEMUX/WholeIsoSeq_Mod.flnc.report.csv

python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq $All_Merged/TOFU/WholeIsoSeq.collapsed.filtered.rep.fa --read_stat $All_Merged/TOFU/WholeIsoSeq.collapsed.read_stat.txt --classify_csv $All_Merged/DEMUX/WholeIsoSeq_Mod.flnc.report.csv --output WholeIsoSeq_Demultiplex_Counts.txt  --primer_names $Demultiplex_Samples

python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq $All_Merged/TOFU/WholeIsoSeq.collapsed.filtered.rep.fa --read_stat $All_Merged/TOFU/WholeIsoSeq.collapsed.read_stat.txt --classify_csv $All_Merged/DEMUX/WholeIsoSeq_Mod.flnc.report.csv --output WholeIsoSeq_Demultiplex_Unlabelled_Counts.txt  -

run_sqanti2 WholeIsoSeq /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/CHECK $All_Merged/TOFU/WholeIsoSeq.collapsed.filtered.gff $All_Merged/DEMUX/WholeIsoSeq_Demultiplex_Counts.txt

############# Alternatively cat multiple refine files and therefore no need to cluster (which takes a lot of time)
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual/REFINE
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/DEMUX
echo "Concatenating:"
ls $REFINE/*flnc.report.csv
# remove all headers before concatenating
awk 'FNR>1' $REFINE/*flnc.report.csv > All.flnc.report.csv
cat $REFINE/K17.flnc.report.csv | head -n 1 > header
cat header All.flnc.report.csv > Final.All.flnc.report.csv
rm header All.flnc.report.csv

Demultiplex_Samples=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Whole_Transcriptome/Demultiplex_IsoSeq_Whole.txt
All_Samples=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples
cd $All_Samples/MOUSE/DEMUX
Rscript /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Whole_Transcriptome/ModifyClassifyReport.R $All_Samples/MOUSE/DEMUX/Final.All.flnc.report.csv $All_Samples/MOUSE/DEMUX/Final_Mod.All.flnc.report.csv

python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq $All_Samples/MOUSE/TOFU/All_Merged.collapsed.filtered.rep.fa --read_stat $All_Samples/MOUSE/TOFU/All_Merged.collapsed.read_stat.txt --classify_csv $All_Samples/MOUSE/DEMUX/Final_Mod1.All.flnc.report.csv --output All_Demultiplex_Counts.txt  --primer_names $Demultiplex_Samples

# run_sqanti2_QC <sample_name> <working_directory> <collapsed.gtf> <count.txt>
run_sqanti2_QC WholeIsoSeq $All_Samples/MOUSE/SQANTI $All_Samples/MOUSE/TOFU/All_Merged.collapsed.filtered.gff $All_Samples/MOUSE/DEMUX/All_Demultiplex_Counts.txt
run_sqanti2_Filter All_Merged $All_Samples/MOUSE/SQANTI
}

#************************************* Rerun with default cupcake parameters (06/01/2021) but simplfiied demutliplexing using custom script
cd $All_Merged; mkdir DEMUX_CUSTOM; cd DEMUX_CUSTOM; mkdir TOFU TAMA_Filter SQANTI
# tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir> <coverage_threshold> <identity_threshold>
# demux <input path read.stat file> <input path of samples file> <path of output>
# run_sqanti2 <sample_name> <working_directory> <collapsed.gtf> <count.txt>
# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
tofu_params WholeIsoSeq $All_Merged/IsoSeq $All_Merged/MAPPING $All_Merged/DEMUX_CUSTOM/TOFU 0.99 0.95
demux $All_Merged/DEMUX_CUSTOM/TOFU/WholeIsoSeq.collapsed.read_stat.txt $SAMPLES_LIST/Demultiplex_IsoSeq_Whole.csv $All_Merged/DEMUX_CUSTOM/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
run_sqanti2 WholeIsoSeq.collapsed.filtered $All_Merged/DEMUX_CUSTOM/SQANTI $All_Merged/DEMUX_CUSTOM/TOFU/WholeIsoSeq.collapsed.filtered.gff $All_Merged/DEMUX_CUSTOM/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
TAMA_remove_fragments $All_Merged/DEMUX_CUSTOM/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite.gtf WholeIsoSeq $All_Merged/DEMUX_CUSTOM/TAMA_Filter

#************************************* Rerun with changed cupcake parameters (06/01/2021)
cd $All_Merged; mkdir DEMUX_CUSTOM_ADJUST; cd DEMUX_CUSTOM_ADJUST; mkdir TOFU TAMA_Filter SQANTI
# tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir> <coverage_threshold> <identity_threshold>
# demux <input path read.stat file> <input path of samples file> <path of output>
# run_sqanti2 <sample_name> <working_directory> <collapsed.gtf> <count.txt>
# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
tofu_params WholeIsoSeq $All_Merged/IsoSeq $All_Merged/MAPPING $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU 0.95 0.95
demux $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.read_stat.txt $SAMPLES_LIST/Demultiplex_IsoSeq_Whole.csv $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
run_sqanti2 WholeIsoSeq.collapsed.filtered $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.filtered.gff $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
TAMA_remove_fragments $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite.gtf WholeIsoSeq $All_Merged/DEMUX_CUSTOM_ADJUST/TAMA_Filter

#************************************* Rerun with changed cupcake parameters to 85% (08/01/2021)
cd $All_Merged; mkdir DEMUX_c85; cd DEMUX_c85; mkdir TOFU TAMA_Filter SQANTI SQANTI_TAMA_FILTER
# tofu <prefix_sample> <input_CLUSTERED_dir> <input_MAPPING_dir> <output_dir> <default/reduced/sigreduced>
# demux <input path read.stat file> <input path of samples file> <path of output>
# run_sqanti2 <sample_name> <working_directory> <collapsed.gtf> <count.txt>
# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
tofu_params WholeIsoSeq $All_Merged/IsoSeq $All_Merged/MAPPING $All_Merged/DEMUX_c85/TOFU sigreduced
demux $All_Merged/DEMUX_c85/TOFU/WholeIsoSeq.collapsed.read_stat.txt $SAMPLES_LIST/Demultiplex_IsoSeq_Whole.csv $All_Merged/DEMUX_c85/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
run_sqanti2 WholeIsoSeq.collapsed.filtered $All_Merged/DEMUX_c85/SQANTI $All_Merged/DEMUX_c85/TOFU/WholeIsoSeq.collapsed.filtered.gff $All_Merged/DEMUX_c85/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
TAMA_remove_fragments $All_Merged/DEMUX_c85/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite.gtf WholeIsoSeq $All_Merged/DEMUX_c85/TAMA_Filter
TAMA_sqanti_filter $All_Merged/DEMUX_c85/TAMA_Filter/WholeIsoSeq.bed $All_Merged/DEMUX_c85/SQANTI WholeIsoSeq.collapsed.filtered_classification.filtered_lite_classification.txt WholeIsoSeq.collapsed.filtered_classification.filtered_lite_classification.gtf WholeIsoSeq.collapsed.filtered_classification.filtered_lite_classification.fasta WholeIsoSeq $All_Merged/DEMUX_c85/SQANTI_TAMA_FILTER

# check number  of isoforms from TAMA filter and retained match that from SQANTI
wc -l $All_Merged/DEMUX_CUSTOM_ADJUST/TAMA_Filter/WholeIsoSeq.bed
wc -l $All_Merged/DEMUX_CUSTOM_ADJUST/TAMA_Filter/WholeIsoSeq_discarded.txt
wc -l $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite_classification.txt

# check demux script does the same thing as Cupcake_Demultiplex.R
# script.py <inputreadstats> <outputfile>
source activate sqanti2_py3
GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
python $GENERALFUNC/Cupcake_Demux_Original.py $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.read_stat.txt $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.Demultiplexed.Cupcake_fake_classify.txt


python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fafq $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.filtered.rep.fa --read_stat $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.read_stat.txt --classify_csv $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.Demultiplexed.Cupcake_fake_classify.txt --output $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.Demultiplexed.Cupcake_Abundance.txt

# Liz's script takes the filtered data, whereas own demux script takes it from the non-filtered data
# the number of isoforms match that of the abundance file (+8 to account for header )
wc -l $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt
wc -l $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.abundance.txt

wc -l $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.Demultiplexed.Cupcake_Abundance.txt
wc -l $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.filtered.abundance.txt

#************************************* Filter SQANTI files using TAMA and feed into IsoformSwitchAnalyze.R
# 1. Classification file with only TAMA selected isoforms
# 2. SQANTI filtered gtf file with only TAMA selected isoforms
# 3. SQANTI filtered fasta file with only TAMA selected isoforms

source activate sqanti2_py3

# Rscript .R <path/tama_filtered_output> <sqanti_filtered_dir> <sqanti_output_name> <output_dir>
Rscript $GENERALFUNC/tama_sqanti_classgtfsubset.R $All_Merged/DEMUX_CUSTOM_ADJUST/TAMA_Filter/WholeIsoSeq.bed $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI WholeIsoSeq.collapsed $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI_TAMA_FILTER
Rscript $GENERALFUNC/tama_sqanti_classgtfsubset.R $All_Merged/DEMUX_CUSTOM/TAMA_Filter/WholeIsoSeq.bed $All_Merged/DEMUX_CUSTOM/SQANTI WholeIsoSeq.collapsed $All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER

# extract fasta sequence based on the pbid (https://www.biostars.org/p/319099/)
# script.py <path/sqanti_filtered.fasta> <path/retained_pbid_tama.txt> <path/output.fasta>
cd $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI_TAMA_FILTER; awk '{ print $4 }' $All_Merged/DEMUX_CUSTOM_ADJUST/TAMA_Filter/WholeIsoSeq.bed | cut -d ";" -f 2  > tama_retained_pbid.txt
python $GENERALFUNC/tama_sqanti_fastasubset.py $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite.fasta $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI_TAMA_FILTER/tama_retained_pbid.txt $All_Merged/DEMUX_CUSTOM_ADJUST/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantifiltered_tamafiltered_classification.fasta

cd $All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER; awk '{ print $4 }' $All_Merged/DEMUX_CUSTOM/TAMA_Filter/WholeIsoSeq.bed | cut -d ";" -f 2  > tama_retained_pbid.txt
python $GENERALFUNC/tama_sqanti_fastasubset.py $All_Merged/DEMUX_CUSTOM/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite.fasta $All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER/tama_retained_pbid.txt $All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantifiltered_tamafiltered_classification.fasta

############## for default setting
mkdir $All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER; cd $All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER
awk '{ print $4 }' $All_Merged/DEMUX_CUSTOM/TAMA_Filter/WholeIsoSeq.bed | cut -d ";" -f 2  > tama_retained_pbid.txt
grep -f tama_retained_pbid.txt $All_Merged/DEMUX_CUSTOM/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite_classification.txt > WholeIsoSeq_sqantifiltered_tamafiltered_classification.txt
head -1 $All_Merged/DEMUX_CUSTOM/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite_classification.txt | cat - WholeIsoSeq_sqantifiltered_tamafiltered_classification.txt > temp && mv temp WholeIsoSeq_sqantifiltered_tamafiltered_headed_classification.txt
grep -f tama_retained_pbid.txt $All_Merged/DEMUX_CUSTOM/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite.gtf > WholeIsoSeq_sqantifiltered_tamafiltered_classification.gtf
source activate sqanti2_py3
python /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/tama_sqanti_fastasubset.py $All_Merged/DEMUX_CUSTOM/SQANTI/WholeIsoSeq.collapsed.filtered_classification.filtered_lite.fasta $All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER/tama_retained_pbid.txt $All_Merged/DEMUX_CUSTOM/SQANTI_TAMA_FILTER/WholeIsoSeq_sqantifiltered_tamafiltered_classification.fasta

######################## Demultiplex only WT samples for differential expression on ages (using 6 samples)
WT_TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/TOFU
WT_SQANTI=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/DEMUX_WT
demux $WT_TOFU/WT_Merged.collapsed.read_stat.txt $SAMPLES_LIST/Demultiplex_IsoSeq_Mouse_WT_Whole.csv $WT_SQANTI/WT.Demultipled_Abundance.txt
run_sqanti2 WT_Merged.collapsed.filtered $WT_SQANTI $WT_TOFU/WT_Merged.collapsed.filtered.gff $WT_SQANTI/WT.Demultipled_Abundance.txt
TAMA_remove_fragments $WT_SQANTI/WT_Merged.collapsed.filtered_classification.filtered_lite.gtf WT_Merged $WT_SQANTI
TAMA_sqanti_filter $WT_SQANTI/WT_Merged.bed $WT_SQANTI WT_Merged.collapsed.filtered_classification.filtered_lite_classification.txt WT_Merged.collapsed.filtered_classification.filtered_lite_classification.gtf WT_Merged.collapsed.filtered_classification.filtered_lite_classification.fasta WT_Merged $WT_SQANTI

################ NCAM1 QC
cd $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU
echo $'PB.16345.39\nPB.16345.63\nPB.16345.57\nPB.16345.73' > NCAM1_identical_subset.txt
python $GENERALFUNC/tama_sqanti_fastasubset.py WholeIsoSeq.collapsed.rep.fasta NCAM1_identical_subset.txt NCAM1_identical.collapsed.rep.fasta
echo $'PB.16345.7\nPB.16345.9' > NCAM1_c99_subset.txt
python $GENERALFUNC/tama_sqanti_fastasubset.py WholeIsoSeq.collapsed.rep.fasta NCAM1_c99_subset.txt NCAM1_c99.collapsed.rep.fasta
echo $'PB.16345.39\nPB.16345.63\nPB.16345.57\nPB.16345.73' > NCAM1_identical_subset.txt
python $GENERALFUNC/tama_sqanti_fastasubset.py WholeIsoSeq.collapsed.rep.fasta NCAM1_identical_subset.txt NCAM1_identical.collapsed.rep.fasta

cd $All_Merged/DEMUX_CUSTOM/TOFU
echo $'PB.16234.33\nPB.16234.54\nPB.16234.50\nPB.16234.64' > NCAM1_identical_subset.txt
python $GENERALFUNC/tama_sqanti_fastasubset.py WholeIsoSeq.collapsed.rep.fasta NCAM1_identical_subset.txt NCAM1_identical.collapsed.rep.fasta


# grep the pbids that are in c95 and not in c99 from cupcake collapse
while read transcript; do
  grep -w $transcript $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/WholeIsoSeq.collapsed.group.txt
done < $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/c95additionaltranscripts.txt > $All_Merged/DEMUX_CUSTOM_ADJUST/TOFU/c95additionaltranscripts_pbid.txt

# check the alignment length and identity of mapped results
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/MAPPING/PAF
htsbox samview -pS All_Merged.clustered.hq.fastq.sam > All_Merged.clustered.hq.fastq.paf
awk -F'\t' '{if ($6="*") {print $0}}' All_Merged.clustered.hq.fastq.paf > All_Merged.allread.clustered.hq.fastq.paf # all reads
awk -F'\t' '{if ($6=="*") {print $0}}' All_Merged.clustered.hq.fastq.paf > All_Merged.notread.clustered.hq.fastq.paf
awk -F'\t' '{if ($6!="*") {print $0}}' All_Merged.clustered.hq.fastq.paf > All_Merged.filtered.clustered.hq.fastq.paf
awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' All_Merged.filtered.clustered.hq.fastq.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > All_Merged_reads_with_alignment_statistics.txt

##################### Checking for hMAPT in WT and TG
MAPT_FUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Isoseq3_QC
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/CLUSTER
SQANTIDEMUX=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/SQANTI_POSTTOFU/
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/MOUSE/TOFU
source $MAPT_FUNC/Mouse_Validation.sh
find_hMAPT $CLUSTER/WT_Merged.clustered.hq.fasta
find_hMAPT $CLUSTER/TG_Merged.clustered.hq.fasta

#Find_Human_Mapt.py <path/input.fasta> <outputname> <path/outputdir>
hMAPT_output_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/QC/Human_Mapt/Sequences
# only in clustered files
python $MAPT_FUNC/Find_Human_Mapt.py $CLUSTER/TG_Merged.clustered.hq.fasta TG_Merged.clustered.hq $hMAPT_output_dir >$hMAPT_output_dir/TG_Merged.clustered.log
python $MAPT_FUNC/Find_Human_Mapt.py $CLUSTER/All_Merged.clustered.hq.fasta All_Merged.clustered.hq $hMAPT_output_dir >$hMAPT_output_dir/All_Merged.clustered.log
python $MAPT_FUNC/Find_Human_Mapt.py $CLUSTER/WT_Merged.clustered.hq.fasta WT_Merged.clustered.hq $hMAPT_output_dir >$hMAPT_output_dir/WT_Merged.clustered.log

# filtered out from cupcake collapse
python $MAPT_FUNC/Find_Human_Mapt.py $SQANTIDEMUX/All_Merged.collapsed.filtered_classification.filtered_lite.fasta All_Merged.sqantifiltered $hMAPT_output_dir > All_Merged.sqantifiltered.log
python $MAPT_FUNC/Find_Human_Mapt.py $TOFU/All_Merged.collapsed.rep.fq All_Merged.collapsed $hMAPT_output_dir > All_Merged.collapsed.log
python $MAPT_FUNC/Find_Human_Mapt.py $TOFU/All_Merged.collapsed.filtered.rep.fa All_Merged.collapsedfiltered $hMAPT_output_dir > All_Merged.collapsedfiltered.log
python $MAPT_FUNC/Find_Human_Mapt.py $TOFU/TG_Merged.collapsed.filtered.rep.fa TG_Merged.collapsedfiltered $hMAPT_output_dir > TG_Merged.collapsedfiltered.log
