################################################################################################
## Define functions for All_Demultiplex_Final_Functions.sh
# 1) run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory>
# 2) run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
# 3) run_REFINE $Sample $Input_LIMA_directory $Output_directory
# 4) run_targeted_REFINE $Input_Pooled_Sample $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
# 5) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# 6) run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
# 7) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
# 8) demux_targeted <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
# 9) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
# 10) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
# 11 run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# 12) run_sqanti2 <input_tofu_prefix> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir>
# 13) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# 14) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
# 15) parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
################################################################################################

#************************************* DEFINE VARIABLES
module load Miniconda2/4.3.21

# Listing versions
source activate isoseq3
ccs --version #ccs 5.0.0 (commit v5.0.0)
lima --version #lima 2.0.0 (commit v2.0.0)
isoseq3 --version #isoseq3 3.4.0 (commit v3.4.0)

source activate sqanti2_py3
samtools --version # echo version
echo "Minimap2 version:" $(minimap2 --version) # echo version
echo "ToFU Cupcake Version"
CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake
head $CUPCAKE/README.md

echo "FASTA SEQUENCE (CLONTECH PRIMERS) FOR NON-MULTIPLEXING"
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
FASTA=$REFERENCE/primer.fasta
cat $FASTA

echo "FASTA SEQUENCE (BARCODED PRIMERS) FOR MULIPLEXING IN TARGETED SEQUENCING"
TARGETED_FASTA=$REFERENCE/targeted.primer.fasta
cat $TARGETED_FASTA

################################################################################################
#*************************************  Isoseq3 [Function 1,2,3,4,5,6]
# run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory>
run_CCS(){
  source activate isoseq3

  cd $3
  echo "Processing Sample $2 from Functions script"
  echo "Processing $1"
  echo "Checking if $2.ccs.bam exists"

    if [ -f $2.ccs.bam ]; then
      echo "$2.ccs.bam file already exists; CCS no need to be processed on Sample $2"
    else
      echo "$2.ccs.bam file does not exist"
      echo "Processing CCS for sample $2"
      # ccs <input.subreads.bam> <output.ccs.bam>
      time ccs $1 $2.ccs.bam --minPasses=1 --min-rq 0.9 --reportFile $2_ccs_report.txt
      echo "CCS for Sample $2 successful"
      ls *$2*
    fi
  source deactivate
}

# run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
# removed --no pbi as this is needed for downstream polishing
run_LIMA(){
  source activate isoseq3

  cd $3
  echo "Processing $1 file for demultiplexing"
    if [ -f $2/$1.fl.json ]; then
      echo "$1.fl.bam file already exists; LIMA no need to be processed on Sample $1"
    elif [ $4 = "multiplex" ]; then
      echo "Multiplex: use Targeted_FASTA"
      #lima <input.ccs.merged.consensusreadset.xml> <input.primerfasta> <output.fl.bam>
      time lima $2/$1.ccs.bam $TARGETED_FASTA $1.fl.bam --isoseq --dump-clips --dump-removed --peek-guess
      echo "lima $1 successful"
      ls $1.fl*
    else
      echo "No-Multiplex: use FASTA"
      time lima $2/$1.ccs.bam $FASTA $1.fl.bam --isoseq --dump-clips --dump-removed
      echo "lima $1 successful"
      ls $1.fl*
    fi
    source deactivate
}

# run_REFINE $Sample $Input_LIMA_directory $Output_directory
run_REFINE(){
  source activate isoseq3

  cd $3
  echo "Processing $1 file for refine"
  if [ -f $1.flnc.bam ]; then
    echo "$1.flnc bam file already exists; Refine no need to be processed on Sample $1"
  else
    #refine --require-polya <input.lima.consensusreadset.xml> <input.primer.fasta> <output.flnc.bam>
    time isoseq3 refine $2/$1.fl.primer_5p--primer_3p.bam $FASTA $1.flnc.bam --require-polya
    echo "refine $1 successful"
    ls $1.flnc*
  fi

  source deactivate
}

# run_targeted_REFINE $Input_Pooled_Sample $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
# Input_config_file = list of samples, with respective barcoded output lima sample name
run_targeted_REFINE(){
  source activate isoseq3

  cd $5
  # extract only the barcoded output lima sample name in config file (note only works if sample has three characters)
  input_sample=$(cat $2 | grep -o "$1.*" | cut -c5-)
  echo "Working with sample $1 = $3.$input_sample"

  time isoseq3 refine $4/$3.$input_sample $TARGETED_FASTA $1.flnc.bam --require-polya
  echo "refine $1 successful"

  source deactivate
}

# merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# aim: merging bam files from refine onwards (similar to run_isoseq3_2_1_merge, but no need to rerun from ccs)
# prerequisite: run ccs, lima, refine for samples
# <input_flnc_bam_dir> = path of flnc bam files generated from refine
# <output_dir> = path of output files from merging
# <output_name> = output name for merged files
# samples.... = list of the sample names
merging_at_refine(){
  module load Miniconda2/4.3.21
  source activate isoseq3
  isoseq3 --version #isoseq3 3.2.2 (commit v3.2.2)

  ###********************* Merging at REFINE
  # Define variable "Merge_Samples" as a list of all samples, in order to find the specified flnc.bam (for dataset create ConsensusReadSet)
  # Define variable "all_flnc_bams" for merged path-directory of all flnc samples (for dataset create ConsensusReadSet)
  Merge_Samples=$(echo "${@:4}")

  echo "Merging flnc of samples $Merge_Samples"
  all_flnc_bams=$(
      for i in ${Merge_Samples[@]}; do
          flnc_bam_name=$(find $1 -name "*.flnc.bam" -exec basename \{} \; | grep ^$i )
          flnc_bam=$(find $1 -name "*.flnc.bam" | grep "$flnc_bam_name" )
          echo $flnc_bam
      done
  )

  cd $2
  printf '%s\n' "${all_flnc_bams[@]}" > $3.flnc.fofn
  cat $3.flnc.fofn

  ###*********************

  isoseq3 cluster $3.flnc.fofn $3.clustered.bam --verbose --use-qvs
  gunzip *.gz*
  # convert clustered.hq.fasta to clustered.hq.fastq for
  source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh
  convert_fa2fq $3.clustered.hq.fasta .

  source deactivate
}

# run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CLUSTER(){
  source activate isoseq3

  cd $3
  echo "Processing $1 file for cluster"
  if [ -f $1.unpolished.bam ]; then
    echo "$1.unpolished.bam file already exists; Cluster no need to be processed on Sample $1"
  else
    # cluster <input.flnc.bam> <output.unpolished.bam>
    time isoseq3 cluster $2/$1.flnc.bam $1.clustered.bam --verbose --use-qvs 2> $1.cluster.log
    echo "cluster $1 successful"
    ls $1.clustered*
  fi

  source deactivate
}


################################################################################################
#************************************* Post_Isoseq3 (Minimap2, Cupcake, Demultiplex) [Function 7,8]
# run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
# Prerequisite: mm10 cage peak
run_map_cupcakecollapse(){

    CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake
    SQANTI2_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    ANNOTATION=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/annotation
    SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence

    source activate sqanti2_py3

    # convert between fasta and fastq for downstream process
    echo "fasta2fastq conversion"
    python $SEQUENCE/fa2fq.py $2/$1.clustered.hq.fasta

    samtools --version # echo version
    echo "Minimap2 version:" $(minimap2 --version) # echo version

    echo "Processing Sample $1 for Minimap2 and sort"
    cd $3 #cd to $MAP directory for output
    minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $REFERENCE/mm10.fa $2/$1.clustered.hq.fastq > $1.sam 2> $1.map.log
    samtools sort -O SAM $1.sam > $1.sorted.sam

    # Alignment stats
    mkdir PAF; cd PAF
    source activate nanopore
    htsbox samview -pS $3/$1.sorted.sam > $1.paf
    awk -F'\t' '{if ($6="*") {print $0}}' $1.paf > $1.allread.paf # all reads
    awk -F'\t' '{if ($6=="*") {print $0}}' $1.paf > $1.notread.paf
    awk -F'\t' '{if ($6!="*") {print $0}}' $1.paf > $1.filtered.paf
    awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $1.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $1"_reads_with_alignment_statistics.txt"
    echo "Number of mapped transcripts to genome:"
    wc -l $1.filtered.paf
    echo "Number of ummapped transcripts to genome:"
    wc -l $1.notread.paf

    echo "Processing Sample $1 for TOFU, with coverage 85% and identity 95%"
    source activate cupcake
    cd $4 #cd to $TOFU directory for output
    collapse_isoforms_by_sam.py -c 0.85 -i 0.95 --input $2/$1.clustered.hq.fastq --fq -s $3/$1.sorted.sam --dun-merge-5-shorter -o $1 &>> $1.collapse.log
    get_abundance_post_collapse.py $1.collapsed $2/$1.clustered.cluster_report.csv &>> $1.abundance.log

    source activate sqanti2_py3
    # convert rep.fq to rep.fa for SQANTI2 input
    seqtk seq -a $1.collapsed.rep.fq > $1.collapsed.rep.fa
    echo "Processing Sample $1 for TOFU successful"
    source deactivate
}


# demux_targeted <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
# read in read.stat file from cupcake collapse output and count abundance of each sample based on CCS ID obtained from each refine report of the sample
demux_targeted(){
    source activate sqanti2_py3
    DEMUXFUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/1_Transcriptome_Annotation/Linux/Targeted_Transcriptome

    # script.R <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
    Rscript $DEMUXFUNCTIONS/Demultiplex_Cupcake.R  $1 $2 $3 $4
}

################################################################################################
#************************************* RNAseq [Function 9, 10]
# run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
# Aim: Individually align RNASeq filtered fastq files to STAR-indexed genome using STAR, and subsequently input SJ.out.bed file into SQANTI2
run_star(){
    source activate sqanti2
    echo "STAR Version"
    STAR --version
    # samtools version: Version: 1.9 (using htslib 1.9)

    # extract official sample name from all_filtered directory
    # extract only files with "fastq.filtered", beginning with sample name, and R1/R2

    # cd to $J20/Tg4510_input_directory
    F_name=$(find $2 -name "*fastq.filtered" -exec basename \{} \; | grep ^$1 | grep "R1" )
    R_name=$(find $2 -name "*fastq.filtered" -exec basename \{} \; | grep ^$1 | grep "R2" )
    # save path directory of files as variable for later mapping
    F_File=$(find $2 -name "$F_name")
    R_File=$(find $2 -name "$R_name")

    # cd to $WKD
    cd $3
    if [ -f $1.SJ.out.bed ]; then
        echo "$1.SJ.out.bedalready exists; STAR Mapping not needed on Sample $1"
    else
        mkdir $1
        cd $1
        echo "Processing Sample $1 for STAR"
        echo "Processing Forward Reads: $F_File"
        echo "Processing Reverse Reads: $R_File"

        STAR_reference_input_directory=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/STAR_main


      # Parameters according to https://github.com/jennylsmith/Isoseq3_workflow/blob/master/shell_scripts/8__STAR_Junctions_ShortReads.sh
      # 2-pass mode alignment with chimeric read detection
      # at least 25 bp of one read of a pair maps to a different loci than the rest of the read pair
      # require 20 pb on either side of chimeric junction
      # include chimeric reads in the output BAM
      # don't include chimeras that map to reference genome contig with Ns
      # --outSAMtype BAM SortedByCoordinate, output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.
      # note STAR indexed with genocde vM22 primary assembly annotation gtf therefore no need to include into command line, otherwise create additional folders
        STAR --runMode alignReads --runThreadN 32 --genomeDir $STAR_reference_input_directory \
        --readFilesIn $F_File $R_File \
        --outSAMtype BAM SortedByCoordinate \
        --chimSegmentMin 25 \
    		--chimJunctionOverhangMin 20 \
    		--chimOutType WithinBAM \
    		--chimFilter banGenomicN \
    		--chimOutJunctionFormat 1 \
    		--twopassMode Basic \
    		--twopass1readsN -1 #use all reads
          #--readFilesCommand zcat \
          #--sjdbGTFfile $4/gencode.vM22.primary_assembly.annotation.gtf \

        # to remove duplicates between samples
        picard MarkDuplicates INPUT=$3/$1".sorted.bam" OUTPUT=$1".sorted.nodup.bam" METRICS_FILE=$1".dup.txt" VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> $1.PicardDup.log

        # rename output files
        mv Aligned.sortedByCoord.out.bam $1.sorted.bam
        mv Log.final.out $1.Log.final.out
        mv Log.out $1.Log.out
        mv Log.progress.out $1.Log.progress.out
        mv SJ.out.tab $1.SJ.out.bed
        mv $1.SJ.out.bed ../
    fi
    source deactivate
}

# mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
# Aim: Merge forward and reverse of all RNASeq filtered fastq, for further downstream alignment to IsoSeq Tofu output using Kallisto
mouse_merge_fastq(){
    R1_READS=($(
        for i in ${SAMPLES_NAMES[@]}; do
            F_name=$(find $1 -name "*fastq.filtered" -exec basename \{} \; | grep ^$i | grep "R1" )
            F_File=$(find $1 -name "$F_name")
            echo $F_File
        done
    ))

    R2_READS=($(
        for i in ${SAMPLES_NAMES[@]}; do
            R_name=$(find $1 -name "*fastq.filtered" -exec basename \{} \; | grep ^$i | grep "R2" )
            R_File=$(find $1 -name "$R_name")
            echo $R_File
        done
    ))

    R1_READS_MERGE=$(echo "${R1_READS[@]}")
    R2_READS_MERGE=$(echo "${R2_READS[@]}")

    echo "Processing R1 READS"
    echo $R1_READS_MERGE
    echo "Processing R2 READS"
    echo $R2_READS_MERGE

    cd $2
    cat $R1_READS_MERGE > $3_R1.fq
    cat $R2_READS_MERGE > $3_R2.fq
}


################################################################################################
#************************************* RNASeq & IsoSeq [Function 11]
# run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# Aim: Align merged RNASeq fastq files to IsoSeq Tofu output fasta using Kallisto
# Prerequisite:
    # run mouse_merge_fastq/any equivalent merging of all RNASeq fastq files (note sample_prefix_output name is the same)
# Output: <output_prefix>.mod.abundance.tsv for input into Sqanti_qc.py (not modified script to take in space delimited rather than tab file)
run_kallisto(){
    source activate sqanti2

    echo "Processing Kallisto for $1"
    cd $4
    kallisto version
    time kallisto index -i $1_Kallisto.idx $2 2> $1_Kallisto.index.log
    time kallisto quant -i $4/$1_Kallisto.idx --fr-stranded $3/$1_R1.fq --rf-stranded $3/$1_R2.fq -o $4 2> $1_Kallisto.quant.log
    mv abundance.tsv $1.abundance.tsv

    # problem: retained co-ordinates, which does not input well into SQANTI2
    echo "Kallisto original $1.abundance.tsv"
    head $1.abundance.tsv
    # solution: retain the PB.ID
    while read line ; do
	    first=$( echo "$line" |cut -d\| -f1 ) # each read line, and remove the first part i.e. PacBio ID
	    rest=$( echo "$line" | cut -d$'\t' -f2-5 ) #save the remaining columns
	    echo $first $rest # concatenate
    done < $4/$1.abundance.tsv > $4/$1.temp.abundance.tsv

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

################################################################################################
#************************************* SQANTI2 [Function 12]
# run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/rnaseq/lncrna>
run_sqanti2(){

    source activate sqanti2_py3

    # variables
    SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    CAGE_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CAGE

    # copy STAR output SJ.out.bed files
    SJ_OUT_BED=($(
        for rnaseq in ${SAMPLES_NAMES[@]}; do
            name=$(find $4 -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
            File=$(find $4 -name "$name")
            echo $File
        done
        ))
    for file in ${SJ_OUT_BED[@]}; do cp $file $7 ;done

    # prepare sqanti
    cd $7
    export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
    python $SQANTI2_DIR/sqanti_qc2.py -v

    # sqanti qc
    echo "Processing Sample $1 for SQANTI2 QC"

    # no kalliso file
    if [ $8 == "rnaseq" ]; then
      python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed --coverage "./*SJ.out.bed" --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt &> $1.sqanti.qc.log

    elif [ $8 == "lncrna" ]; then
      echo "Processing with lncRNA.gtf for genome annotation "
      python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM25.long_noncoding_RNAs.gtf $REFERENCE/mm10.fa --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed --coverage "./*SJ.out.bed" --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --skipORF --fl_count $6 &> $1.sqanti.qc.log

    elif [ $8 == "genome" ]; then
      echo "Processing with gencode.vM22.annotation.gtf for genome annotation "
      python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed --coverage "./*SJ.out.bed" --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --expression $5 --fl_count $6 &> $1.sqanti.qc.log
    else
      echo "8th argument required"
    fi

    echo "Processing Sample $1 for SQANTI2 filter"
    python $SQANTI2_DIR/sqanti_filter2.py $1"_classification.txt" $1"_corrected.fasta" $1"_corrected.gtf" -a 0.6 -c 3 &> $1.sqanti.filter.log

    # remove temp SJ.out bed files
    rm *SJ.out.bed
    source deactivate
}



################################################################################################
#************************************* TAMA [Function 13,14]
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

# TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <sqanti_output_junc> <output_prefix_name> <output_dir>
TAMA_sqanti_filter(){
  source activate sqanti2_py3
  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
# Rscript .R <path/tama_filtered_output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_junc_txt> <output_prefix_name> <output_dir>
  Rscript $GENERALFUNC/TAMA/tama_sqanti_classgtfsubset.R $1 $2 $3 $4 $6 $7 $8
  # extract fasta sequence based on the pbid (https://www.biostars.org/p/319099/)
  # script.py <path/sqanti_filtered.fasta> <path/retained_pbid_tama.txt> <path/output.fasta>
  cd $8
  awk '{ print $4 }' $1| cut -d ";" -f 2  > tama_retained_pbid.txt
  python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $2/$5 $8/tama_retained_pbid.txt $8/$7"_sqantifiltered_tamafiltered_classification.fasta"

  SQ_Report=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2/utilities/SQANTI_report2.R
  Rscript $SQ_Report $7"_sqantitamafiltered.classification.txt" $7"_sqantitamafiltered.junction.txt" $>> $7"_sqantitamafiltered.sq_qc.log"
}


################################################################################################
#************************************* QC [Function 15]
# parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
parse_stats_per_sample(){

    # variable
    CCS_dir=$1
    LIMA_dir=$2
    output_name=$3

    FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
    source activate sqanti2_py3
    python $FUNCTIONS/IsoSeq_QC/CCS.py $CCS_dir "" $3
    python $FUNCTIONS/IsoSeq_QC/LIMA.py $LIMA_dir "" $3
    source deactivate
}
