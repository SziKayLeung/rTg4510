################################################################################################
# 1) run_QC <sample> <sequencing_summary> <bam_input> <output_dir>


################################################################################################


module load Miniconda2/4.3.21

# 1) run_merge <raw_directory> <sample_output_name>
# output: <sample_output_name>.merged.fastq in $RAWDATA (defined in functions script)
run_merge(){

  source activate nanopore
  
  echo "Merging following fastq files"
  FASTQ=$(ls $1/*) 
  echo ${FASTQ}

  cat ${FASTQ} > $2
  echo "Merge of Samples successful: output to $2"
  
}

merge_fastq_across_samples(){
  # variables 
  gval=$1
  input_dir=$2
  output_dir=$3
  
  echo "Merging ${gval}"
  
  fastq=$(for f in ${input_dir}/*; do ls $f/*${gval}*; done)
  cat ${fastq} > ${output_dir}/${gval}_merged.fastq
}

# 2) run_QC <sample> <sequencing_summary> <bam_input> <output_dir>
run_QC(){
    # variables
    sample=$1
    sequencing_summary=$2
    bam_input=$3
    output_dir=$4

    source activate nanopore
    echo "Processing: $1"
    cd $output_dir
    pycoQC --summary_file $sequencing_summary --bam_file $bam_input -o $sample"_QC.html"
    Rscript ${MINIONQC} -i $sequencing_summary -s TRUE -o $output_dir
    source deactivate
}

# 3) run_porechop <raw.fastq.gz> <output_dir>
# input: <raw>.fastq.gz
# output: <output_directory> with barcodes demultiplexed
run_porechop(){

    source activate nanopore
    sample=$(basename $1)
    
    echo "Processing Sample $sample for Porechop"
    python ${PORECHOP} -i $1 -b $2 --format fastq --threads 16 \
      --check_reads 1000 \
      --discard_middle \
      --end_size 100 \
      --min_trim_size 15 \
      --extra_end_trim 1 \
      --end_threshold 75 \
      --verbosity 2
}


# 4) post_porechop_run_cutadapt <input_fastq> <output_dir>
post_porechop_run_cutadapt(){
  
  input_dir=$(dirname $1)
  name=$(basename $1 .fastq)
  
  source activate nanopore 
  
  # requires fasta files for downstream
  echo "Converting $1 to fasta"
  seqtk seq -a $1 > ${input_dir}/${name}.fasta
  
  # subset fasta file to polyA and polyT fasta (i.e. reads ending with PolyA and starting with polyT)
  # reads that end with AAAAAAAAAA = plus reads 
  # reads that start with TTTTTTTTTT = minus reads (need to be reverse complemented)
  echo "Subsetting fasta to polyA and polyT sequences"
  python ${SUBSETPOLYTAILS} --fa ${input_dir}/${name}.fasta --o_name ${name} --o_dir $2
  
  # working in output directory
  cd $2
  
  # reverse complement minus reads (reads ending with polyT)
  seqtk seq -r ${name}_PolyT.fasta > ${name}_PolyT_rev.fasta
  
  # use cutadapt package to trim polyA
  echo "Remove polyA sequences using cutadapt"
  cutadapt -a "A{60}" ${name}_PolyA.fasta -o ${name}_PolyA_cutadapted.fasta &> ${name}_polyA_cutadapt.log
  cutadapt -a "A{60}" ${name}_PolyT_rev.fasta -o ${name}_PolyT_rev_cuptadapted.fasta &> ${name}_polyT_cutadapt.log
  
  # concatenated reverse minus polyT and polyA reads
  cat ${name}_PolyA_cutadapted.fasta ${name}_PolyT_rev_cuptadapted.fasta > ${name}_combined.fasta
  
  source deactivate
}


# 6) run_minimap2 <input_fasta> <output_dir>
# Aim: Align reads from trimming, filtering to genome of interest using Minimap2
# Input: <sample_name>_combined_reads.fasta
# Output: <sample_name>_combined_reads.sam, <sample_name>_Minimap2.log
run_minimap2(){

	source activate nanopore
	
	name=$(basename $1 .fasta)
	echo "Aligning ${name} using Minimap2"
	
	minimap2 -t 46 -ax splice ${GENOME_FASTA} $1 > $2/${name}.sam 2> $2/${name}_minimap2.log
  samtools sort -O SAM $2/${name}.sam > $2/${name}_sorted.sam

  source deactivate
}


# run_transcriptclean <input_sam> <output_dir>
run_transcriptclean(){
  
  source activate sqanti2_py3
   
  name=$(basename $1 _merged_combined_sorted.sam)
  echo "TranscriptClean ${name}"
  
  cd $2; mkdir -p ${name}
  cd $2/${name}
  python ${TCLEAN} --sam $1 --genome ${GENOME_FASTA} --outprefix $2/${name}/${name} --tmpDir $2/${name}/${name}_tmp
}

# filter_alignment <input_name> <input_mapped_dir>
filter_alignment(){
  
  source activate nanopore
  
  cd $2
  echo "Converting bam to sam and sort"
  samtools view -h $1.bam > $1.sam
  samtools bam2fq $1.bam| seqtk seq -A > $1.fa
  samtools sort -O SAM $1.sam > $1.sorted.sam

  # Alignment stats
  # Use the inforation in the paf file to create a new file where the columns correspond to the following: 
    #col1: name of the nanopore read 
    #col2: name of the sequence where nanopore read aligns (target sequence)
    #col3: start position of the alignment on the target sequence 
    #col4: length of the original nanopore read 
    #col5: length of the aligned part of the nanopore read  
    #col6: fraction of the aligned part of the nanopore read over the orginal length 
    #col7: fraction of the aligned part of the target sequence over the orginal length of the target sequence
    #col8: strand where the nanopore read aligns
    #col8: number of matched nucleotides of the nanopore read alignment on the target sequence
    #col9: identity (percentage of matched nucleotides over the aligned length of the nanopore read)
    #col10: number of mismatches of the nanopore read alignment on the target sequence
    #col11: number of insertions of the nanopore read alignment on the target sequence
    #col12: number of deletions of the nanopore read alignment on the target sequence
  
  echo "Dissecting alignment statistics"
  mkdir -p PAF; cd PAF
  htsbox samview -pS $2/$1.sorted.sam > $1.paf
  awk -F'\t' '{if ($6!="*") {print $0}}' $1.paf > $1.filtered.paf
  awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $1.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $1"_mappedstats.txt"
  ## filter based on alignable length (>0.85) and identity (>0.95)
  awk -F'\t' '{if ($6>=0.85 && $8>=0.95) {print $1}}' $1"_mappedstats.txt" > $1_filteredreads.txt

  source activate sqanti2
  picard FilterSamReads I=$2/$1.bam O=$2/$1.filtered.bam READ_LIST_FILE=$2/PAF/$1_filteredreads.txt FILTER=includeReadList &> $2/PAF/$1.picard.log
  
  source activate nanopore
  samtools bam2fq $2/$1.filtered.bam| seqtk seq -A > $2/$1.filtered.fa
  samtools sort -O SAM $2/$1.filtered.bam -o $2/$1.filtered.sorted.bam
  
  # https://bioinformatics.stackexchange.com/questions/3380/how-to-subset-a-bam-by-a-list-of-qnames
  #source activate nanopore
  #samtools view $2/$1.bam | grep -f $1_filteredreads.txt > $1.filtered.sam
  #samtools view -bS $1.filtered.sam > $1.filtered.bam
  #samtools bam2fq $2/$1.filtered.bam| seqtk seq -A > $2/$1.filtered.fa

}


# run_talon_label <input_sam> <output_dir>
run_talon_label(){
  
  name=$(basename $1 _clean.sam)
  echo "Label ${name} for TALON"
  
  source activate sqanti2_py3
  talon_label_reads --f $1 --g ${GENOME_FASTA} --t 1 --ar 20 --tmpDir=$2/${name}_label_reads --o $2/${name}
  
  # convert sam to fasta
  source activate nanopore
  cd $2 
  samtools view -bS ${name}_labeled.sam > ${name}_labelled.bam
  samtools bam2fq ${name}_labelled.bam | seqtk seq -A > ${name}_labelled.fasta
  
  source deactivate 
}


# run_talon <config_file> 
run_talon(){
  
  source activate sqanti2_py3
  
  output_dir=$(dirname $1)
  
  cd $output_dir
  echo "Running talon using ${TALON_DB}/${SPECIES}_talon.db"  
  talon --f $1 --db ${TALON_DB}/${SPECIES}_talon.db --build ${SPECIES} --o ${NAME}
}

#post_talon <output_dir>
post_talon(){
  
  source activate sqanti2_py3
  cd $1; mkdir -p 1_unfiltered 2_filtered

  talon_cmd=$(echo --db ${TALON_DB}/${SPECIES}_talon.db -a ${SPECIES}_annot --build ${SPECIES})
  
  # Unfiltered dataset 
  echo "Generate abundance and gtf file for unfiltered dataset"
  cd $1/1_unfiltered
  talon_abundance ${talon_cmd} --o ${NAME}_unfiltered
  talon_create_GTF ${talon_cmd} --o ${NAME}_unfiltered
  
  echo "Filtering transcripts using all samples"
  cd $1/2_filtered
  datasets=$(echo ${TALON_FILTER_SAMPLES[@]} | tr ' ' ,)
  talon_filter_transcripts --db ${TALON_DB}/${SPECIES}_talon.db --datasets ${datasets} -a ${SPECIES}_annot --maxFracA 0.5 --minCount 5 --minDatasets 2 --o $1/2_filtered/filtered_transcripts.csv
  
  # Filtered dataset 
  echo "Generate abundance and gtf file for filtered dataset"
  talon_abundance ${talon_cmd} --whitelist $1/2_filtered/filtered_transcripts.csv --o ${NAME}_filtered
  talon_create_GTF ${talon_cmd} --whitelist $1/2_filtered/filtered_transcripts.csv --o ${NAME}_filtered
  
  source deactivate
}


# run_sqanti3 <gtf> <output_dir>
run_sqanti3(){
  
  source activate sqanti2_py3
  
  name=$(basename $1 .gtf)

  cd $2
  
  # sqanti qc
  echo "Processing Sample ${name} for SQANTI3 QC"
  python $SQANTI3_DIR/sqanti3_qc.py -v
  echo ${GENOME_GTF}
  echo ${GENOME_FASTA}
  
  python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $1 ${GENOME_GTF} ${GENOME_FASTA} \
  --cage_peak ${CAGE_PEAK} \
  --polyA_motif_list ${POLYA} \
  --genename --isoAnnotLite --gff3 $GFF3 --report pdf &> ${name}.sqanti.qc.log
  
  echo "Processing Sample ${name} for SQANTI filter"
  python $SQANTI3_DIR/sqanti3_RulesFilter.py ${name}"_classification.txt" ${name}"_corrected.fasta" ${name}"_corrected.gtf" -a 0.6 -c 3 &> ${name}.sqanti.filter.log
  
  Rscript ${SQ_Report} $2/${name}"_classification.txt" $2/${name}"_junctions.txt"
  Rscript ${SQ_Report} $2/${name}"_classification.filtered_lite_classification.txt" $2/${name}"_classification.filtered_lite_junctions.txt"
  
  source deactivate
}
