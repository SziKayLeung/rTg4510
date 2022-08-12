# Szi Kay Leung (sl693@exeter.ac.uk)
################################################################################################
## Define functions for <Sample>_Nanopore.sh
# 1) run_merge <input_fasta_dir> <sample_output_name>
# 2) run_nanofilt <sample_output_name_merged.fastq> <sample_output_name> <root_dir>
# 3) run_porechop <sample> <root_directory> <type>
# 4) post_porechop <sample_output_name> <input_directory>
# 5) run_cutadapt_and_combine <sample_output_name> <input/output_directory>
# 6) run_minimap2 <sample_name> <input_directory> <output_directory>
# 7) filter_mapped_reads <sample_name> <input/output_directory>
# 8) run_tama_collapse <input_drectory/sample_name> <sample_output_tama_collapse_file_prefix> <output_directory> <type>
# 9) Modify_genomegtf_TAMAinput <genocode.gtf path> <list_of_gene_transcript_id_path> <output_name> <output_dir>
  # Prequisite: run "list_mm10_gene_trascript.R" to generate <list_of_gene_transcript_id_path>
# 10) run_tama_merge <genome_bed12> <sample_output_tama_collapse_bed> <sample_output_tama_merge_file_prefix> <output_directory>

################################################################################################
# 13/06/2019: Define Functions for Nanpore Analysis
# 14/04/2020: Estabilshed pipeline from WTAC

#************************************* DEFINE FUNCTIONS
module load Miniconda2

# 1) run_merge <input_fasta_dir> <sample_output_name>
# output: <sample_output_name>.merged.fastq in $RAWDATA (defined in functions script)
run_merge(){

  source activate nanopore
  
  echo "Merging following fastq files in $1"
  cd $1
  zcat *.gz > $RAW_ROOT_DIR/$2"_Merged.fq"
  source deactivate

}

# 2) run_nanofilt <sample_output_name_merged.fastq> <sample_output_name> <root_dir>
# output: <sample_output_name>_filtered_pass_reads.fastq, <sample_output_name>_NanoFilt.log
run_nanofilt(){

    source activate nanopore
    mkdir -p $3/1_nanofilt; cd $3/1_nanofilt

    echo "Processing: $1"
    cat $1 | NanoFilt -q 7 > $/$2_filtered_pass_reads.fastq
    mv NanoFilt.log $3/$2_NanoFilt.log

    source deactivate
    
}

# 2) run_QC <sample> <sequencing_summary> <bam_input> <output_dir>
run_QC(){
    # variables
    sample=$1
    sequencing_summary=$2
    bam_input=$3
    output_dir=$4
    minionQC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/minion_qc/MinIONQC.R

    source activate nanopore
    echo "Processing: $1"
    cd $output_dir
    pycoQC --summary_file $sequencing_summary --bam_file $bam_input -o $sample"_QC.html"
    Rscript $minionQC -i $sequencing_summary -s TRUE -o $output_dir
    source deactivate
}


# on_target_rate <sample> <input_basecalled_dir> <input_mapped_dir> <input_bed_file> <output_directory>
on_target_rate(){
    # variables
    sample=$1
    input_basecalled_dir=$2
    input_mapped_dir=$3
    input_bed_file=$4
    output_dir=$5

    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake

    source activate sqanti2_py3
    cd $output_dir

    export PYTHONPATH=$PYTHONPATH:$SEQUENCE

    echo "Processing combined.fasta: $input_basecalled_dir/$sample"_combined_reads.fasta""
    echo "Processing input_mapped.fasta.sam: $input_mapped_dir/$sample"_combined_sorted_reads.sam""
    echo "Writing output: $sample"
    python $CUPCAKE/targeted/calc_probe_hit_from_sam.py $input_bed_file $input_basecalled_dir/$sample"_combined_reads.fasta" $input_mapped_dir/$sample"_combined_sorted_reads.sam" --start_base 0 --end_base 0 --o $sample.fasta.sam.probe_hit.txt
    source deactivate
}


# 3) run_porechop <sample> <root_directory> <type>
# input: <sample_output_name>_filtered_pass_reads.fastq
# output: <sample_output_name>_filtered_reads_choped.fq, <sample_output_name>_output_adaptor_alignment_stat, <sample_output_name>_grep_output_adaptor_alignment_stat
run_porechop(){
  
    mkdir -p $2/2_demultiplex; cd $2/2_demultiplex
    mkdir -p $2/2_demultiplex/stats
    mkdir -p $2/2_demultiplex/chopped

    source activate nanopore
    
    echo "Processing Sample $1 for Porechop"
    
    if [ $3 == "Whole" ]; then
      python $PORECHOP_WHOLE -i $2/1_nanofilt/$1"_filtered_pass_reads.fastq" -o $2/2_demultiplex/$1"_filtered_reads_choped.fq" --format fastq --verbosity 4 --threads 16 --check_reads 1000 --discard_middle --end_size 100 --min_trim_size 15 --extra_end_trim 1 --end_threshold 75 > $3/$2"_output_adaptor_alignment_stats"
  
    elif [ $3 == "Targeted" ]; then
      # the barocode option in Porechop difficult to use as not able to diffrentiate between the plus and minus strand 
      # also issues and unclear instructions --> need forward and reverse in the name when defining the barcodes    
      python $PORECHOP_TARGETED -i $2/1_nanofilt/$1"_filtered_pass_reads.fastq" -o $2/2_demultiplex/$1"_filtered_reads_choped.fq" --format fastq --verbosity 4 --threads 16 --check_reads 1000 --discard_middle --end_size 100 --min_trim_size 15 --extra_end_trim 1 --end_threshold 75 > $3/$2"_output_adaptor_alignment_stats"
  
    else
      echo "4th argument required as Whole or Targeted"
    fi
  
    #echo "Processing Sample $2 for Porechop Arrangment (WTAC)"
    #sed -i 's/adaptors_at_the_end/adaptors_at_the_end: /g' $2_output_adaptor_alignment_stats
    grep -e "^read_name:" -e "^adaptors_at_the_start:" -e "^adaptors_at_the_end:" $2/2_demultiplex/$1"_output_adaptor_alignment_stats" > $2/2_demultiplex/$1"_grep_output_adaptor_alignment_stats"
    
    mv *filtered_reads_choped* chopped
    mv *alignment_stats* stats
  
    source deactivate

}

# 4) post_targeted_porechop <sample> <root_directory>
post_targeted_porechop(){

    source activate nanopore

  	cd $2/2_demultiplex/stats
  	# Note this is not the original WTAC script, but simplified and allows user input
  	# parse_the_output_adaptor_alignment_stats_file <sample_name>
    echo "Parsing through $1_grep_output_adaptor_alignment_stats: Output as $1_flatten_grep_output_adaptor_alignment_stats"
  	python $PARSEPORECHOP $2/$1_grep_output_adaptor_alignment_stats $2/$1_flatten_grep_output_adaptor_alignment_stats 
   
    #Rscript TargetedTranscriptome_SampleDemultiplex.R <porechop_flatten_stats> <output_dir>
    mkdir -p $2/2_demultiplex/$1"_Demultiplex"; cd $2/2_demultiplex/$1"_Demultiplex"
    Rscript $ONTDEMUX $2/$1_flatten_grep_output_adaptor_alignment_stats  $2/2_demultiplex/$1"_Demultiplex" &> $1"_demultiplex.log"
   
    #Extract sequence from porechop.trimmed.fastq using the names (identifiers)
    cd $2/2_demultiplex/$1"_Demultiplex"
    for bc in $(seq 1 $MAX_BARCODE); do
         
      echo "Processing BC"$bc"_Plus_Strand_Reads.txt to BC"$bc"_Plus_Strand_Reads.fq" 
      seqtk subseq $2/2_demultiplex/chopped/$1_filtered_reads_choped.fq "BC"$bc"_Plus_Strand_Reads.txt" > "BC"$bc"_Plus_Strand_Reads.fq"
      
      echo "Processing BC"$bc"_Minus_Strand_Reads.txt to BC"$bc"_Minus_Strand_Reads.fq" 
      seqtk subseq $2/2_demultiplex/chopped/$1_filtered_reads_choped.fq "BC"$bc"_Minus_Strand_Reads.txt" > "BC"$bc"_Minus_Strand_Reads.fq"
      
      echo "Reverse complement BC"$bc"_Minus_Strand_Reads.fq to BC"$bc"_Reverse_Minus_Strand_Reads.fq" 
      seqtk seq -r "BC"$bc"_Minus_Strand_Reads.fq" > "BC"$bc"_Reverse_Minus_Strand_Reads.fq"
    done
    
   	# All Reads by strand
   	cd $2/2_demultiplex
    seqtk subseq $2/2_demultiplex/chopped/$1"_filtered_reads_choped.fq" $2/2_demultiplex/$1"_Demultiplex"/All_Plus_Reads.txt > $1"_plus_reads.fastq"
  	seqtk subseq $2/2_demultiplex/chopped/$1"_filtered_reads_choped.fq" $2/2_demultiplex/$1"_Demultiplex"/All_Minus_Reads.txt > $1"_minus_reads.fastq"
    seqtk seq -r $1"_minus_reads.fastq" > $1"_reverse_minus_reads.fastq"
    
}



# 5) run_cutadapt_and_combine <sample> <root_directory> <demux/nondemux>
# Aim: trim 40As at the end of reverse complement minus reads and plus reads
# Input: <sample_output_name>reverse_minus_reads.fastq; <sample_output_name>_plus_reads_fastq
# Output: <sample_output_name>_combined_reads.fastq/fasta
run_cutadapt_and_combine(){

  source activate nanopore
  
  mkdir -p $2/3_cutadapt_merge; cd $2/3_cutadapt_merge

  if [ $3 == "demux" ]; then
    reverse_minus=$2/2_demultiplex/$1"_Reverse_Minus_Strand_Reads.fq"
    plus_reads=$2/2_demultiplex/$1"_Plus_Strand_Reads.fq"
  elif [ $3 == "nondemux" ]; then 
    reverse_minus=$2/2_demultiplex/$1"_reverse_minus_reads.fastq"
    plus_reads=$2/2_demultiplex/$1"_plus_reads.fastq"
  else
    echo "demux or nondemux required for 3rd parameter"
  fi
  
  echo "Using:" 
  ls $2/2_demultiplex/$reverse_minus
  ls $2/2_demultiplex/$plus_reads

	# use cutadapt package to trim polyA
	cutadapt -a "A{60}" -o $1_nopolya_reverse_minus_reads.fastq $reverse_minus &> $1_minus_cutadapt.log
	cutadapt -a "A{60}" -o $1_nopolya_plus_reads.fastq $plus_reads &> $1_plus_cutadapt.log

	# concatenated reverse minus and positive reads
	cat $1_nopolya_reverse_minus_reads.fastq $1_nopolya_plus_reads.fastq > $1_combined_reads.fastq

	# convert concatenated fastq into fasta reads
	seqtk seq -a $1_combined_reads.fastq > $1_combined_reads.fasta

  source deactivate

}

# 6) run_minimap2 <sample_name> <root_directory> <batch>
# Aim: Align reads from trimming, filtering to genome of interest using Minimap2
# Input: <sample_name>_combined_reads.fasta
# Output: <sample_name>_combined_reads.sam, <sample_name>_Minimap2.log
run_minimap2(){

	source activate nanopore
	
	mkdir -p $2/4_minimap $2/4_minimap/$3
	cd $3/4_minimap/$3
	
	echo "Alignment using Minimap2: $1"
	
	minimap2 -t 46 -ax splice $GENOME_FASTA $2/2_demultiplex/$3"_Demultiplex"/$1"_combined_reads.fasta" > $1"_combined_reads.sam" 2> $1"_Minimap2.log"
  samtools sort -O SAM $1"_combined_reads.sam" > $1"_combined_sorted_reads.sam"

  source deactivate
  
}

# run_transcriptclean <sample> <root_dir> <batch>
run_transcriptclean(){
  
  mkdir -p $2/5_tclean $2/5_tclean/$3
	cd $2/5_tclean/$3
  
  source activate sqanti2_py3
  python $TRANSCRIPTCLEAN --sam $2/4_minimap/$3/$1"_combined_reads.sam" --genome $GENOME_FASTA --outprefix $2/5_tclean/$3/$1 --tmpDir $2/5_tclean/$3/$1"_transcriptcleantmp"
  
}

# create_talon_db <talon_database_name>
create_talon_db(){
  
  mkdir -p $2/6_talon; cd $2/6_talon
  
  source activate sqanti2_py3
  talon_initialize_database --f $GENOME_GTF --a $1"_annot" --g $1 --o $1"_talon" &> talon_init.log
  
}

# run_talon_label <sample> <root_dir> <batch>
run_talon_label(){
  
  mkdir -p $2/6_talon/1_label $2/6_talon/1_label/$3
  
  source activate sqanti2_py3
  
  # label reads
  talon_label_reads --f $2/5_tclean/$3/$1"_clean.sam" --g $GENOME_FASTA --t 16 ar --ar 20 --o $2/6_talon/1_label/$3/$1 --deleteTmp
  
  # convert sam file to fasta
  cd $2/6_talon/1_label/$3
  samtools view -bS $1"_labeled.sam" > $1"_labelled.bam"
  samtools bam2fq $1"_labelled.bam" | seqtk seq -A > $1"_labelled.fasta"
  
}

# run_talon_and_quantify <name> <root_dir>
run_talon_and_quantify(){
  
  source activate sqanti2_py3
  # run talon
  talon --f $TALON_CONFIG --db $TALON_NAME"_talon.db" --build $TALON_NAME --o $1
  
  # generate abundance 
  talon_abundance --db $TALON_NAME"_talon.db" -a $TALON_NAME"_annot" --build $TALON_NAME --o $1"_unfiltered"
  
  # create gtf 
  talon_create_GTF --db $TALON_NAME"_talon.db" -a $TALON_NAME"_annot" --build $TALON_NAME --o $1"_unfiltered"
  
  mkdir -p $2/6_talon/2_talon_full
  mv *$1"_unfiltered"* $2/6_talon/2_talon_full
  
}


# talon_filter_quantify <name> <root_dir>
talon_filter_quantify(){
  
  mkdir -p $2/6_talon/3_talon_filter
  
  source activate sqanti2_py3
  
  # talon filtering 
  talon_filter_transcripts \
    --db $2/6_talon/mm10_talon.db \
    --datasets $TALON_FILTER_SAMPLES \
    -a mm10_annot --maxFracA $TALON_MAX_A_PARAM --minCount $TALON_MIN_COUNT_PARAM --minDatasets $TALON_MIN_DATASET_PARAM \
    --o $2/6_talon/3_talon_filter/filtered_transcripts.csv
  
  # generate abundance
  wlist=$2/6_talon/3_talon_filter/filtered_transcripts.csv
  talon_abundance --db $TALON_NAME"_talon.db" --whitelist $wlist -a $TALON_NAME"_annot" --build $TALON_NAME --o $1"_filtered"
  
  # generate gtf
  talon_create_GTF --db $TALON_NAME"_talon.db" --whitelist $wlist -a $TALON_NAME"_annot" --build $TALON_NAME --o $1"_filtered"
  
  mv *$1"_filtered"* $2/6_talon/3_talon_filter

}


################################################################################################
#************************************* SQANTI3 [Function 12]

# run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna> <input_dir> <dataset>
run_sqanti3(){

    source activate sqanti2_py3
    
    # variable 
    sample=$1
    gtf=$3/$sample.gtf
    #abundance=
    
    # create directory
    mkdir -p $3/7_sqanti3 
    mkdir -p $3/7_sqanti3/$4 $3/7_sqanti3/$4/$2 
    cd $3/7_sqanti3/$4/$2
    
    # copy STAR output SJ.out.bed files
    SJ_OUT_BED=($(
        for rnaseq in ${RNASEQ_SAMPLES_NAMES[@]}; do
            name=$(find $RNASEQ_MAPPED_DIR -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
            File=$(find $RNASEQ_MAPPED_DIR -name "$name")
            echo $File
        done
        ))
    for file in ${SJ_OUT_BED[@]}; do cp $file $3/7_sqanti3/$4/$2;done
  
    
    # sqanti qc
    echo "Processing Sample $sample for SQANTI3 QC"
    python $SQANTI3_DIR/sqanti3_qc.py -v
    echo $GENOME_GTF
    echo $GENOME_FASTA
    
    if [ $2 == "basic" ]; then
      echo "Processing basic commands"
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $gtf $GENOME_GTF $GENOME_FASTA --cage_peak $CAGE_PEAK --polyA_motif_list $POLYA \
      --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report pdf &> $sample.sqanti.qc.log
    else
      echo "2nd argument required"
    fi

    echo "Processing Sample $sample for SQANTI2 filter"
    python $SQANTI3_DIR/sqanti3_RulesFilter.py $sample"_classification.txt" $sample"_corrected.fasta" $sample"_corrected.gtf" -a 0.6 -c 3 &> $1.sqanti.filter.log
    
    if [ $2 != "basic" ]; then 
      # remove temp SJ.out bed files
      rm *SJ.out.bed
    fi
    
    source deactivate
}



# run_cpat <input_fasta> <species="Mouse"/"Human"> <output_name> <output_dir>  
run_cpat(){
  source activate sqanti2_py3
  
  #variables 
  input_fasta=$1
  species=$2
  output_name=$3
  output_dir=$4
  
  cd $output_dir
  cpat.py -x $REFERENCE/CPAT/$species"_Hexamer.tsv" -d $REFERENCE/CPAT/$species"_logitModel.RData" \
  -g $input_fasta --min-orf=50 --top-orf=50 -o $output_name 2> $output_name.e

}

TAMA_remove_fragments(){

    TAMAFUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/TAMA
    TAMA_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama/tama_go/filter_transcript_models
    source activate sqanti2_py3
    cd $3

    # convert gtf to bed12
	  gtfToGenePred $1 $2.genepred
	  genePredToBed $2.genepred > $2.bed12
    awk -F'\t' '{print $1,$2,$3,$4,"40",$6,$7,$8,"255,0,0",$10,$11,$12}' $2.bed12| sed s/" "/"\t"/g|sed s/",\t"/"\t"/g|sed s/",$"/""/g > Tama_$2.bed12

    # Rscript script.R <name>.bed12 <input_dir>
	  Rscript $TAMAFUNCTIONS/TAMA_Merge_Prepare.R Tama_$2 $3
    source activate sqanti2
	  python $TAMA_DIR/tama_remove_fragment_models.py -f Tama_$2_mod.bed12 -o $2 -m $4

	  rm *Tama*

    echo "Number of isoforms filtered by TAMA:"
    wc -l $2"_discarded.txt"

	  source deactivate
}


convert_gtf_bed12(){

  	source activate nanopore
    cd $2
    sample=$(basename "$1" | cut -d "." -f 1 )
    echo $sample
    gtfToGenePred $sample.gtf $sample.genePred
    genePredToBed $sample.genePred $sample.bed12
    sort -k1,1 -k2,2n $sample.bed12 > $sample"_sorted.bed12"
    rm $sample.genePred $sample.bed12
    # Rscript script.R <input.classfile> <input_output_dir of bed file> <prefix>
    #Rscript $GENERAL/annotate_uscs_tracks.R $SQANTI $WKD $sample 
    #bedToBigBed -extraIndex=name $sample $sample"_Modified.bed12" $REFERENCE/mm10.chrom.filtered.sizes $sample.bb
    
    #source deactivate
    #./bedToBigBed -as=bedExample2.as -type=bed9+3 -extraIndex=name $sample"_Modified.bed12" $REFERENCE/mm10.chrom.filtered.sizes $sample.bb
}