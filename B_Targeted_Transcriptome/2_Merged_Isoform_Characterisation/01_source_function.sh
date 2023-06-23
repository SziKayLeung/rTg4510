# source functions
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
FICLE_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/FICLE
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
export PATH=$PATH:${LOGEN_ROOT}/merge_characterise_dataset
export PATH=$PATH:${FICLE_ROOT}
export PATH=$PATH:${FICLE_ROOT}/reference


#run_gff_compare <output_name> <ref_gtf> <2nd_gtf> <root_dir>
run_gff_compare(){
  
  mkdir -p $4/1_gffcompare; cd $4/1_gffcompare
  cp $2 . 
  cp $3 . 
  
  source activate sqanti2_py3 
  gffcompare -r $2 $3 -o $1
  
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
  
 
# identify_common_transcripts <root_dir> <targetgenelist> <samplelist>
# Custom script to identify_common_targeted_transcripts for merged output
identify_common_transcripts(){
  
  mkdir -p $1/2_common_transcripts
  
  source activate sqanti2_py3
  python $COMMONSC \
  --iso $ISOSEQ_SQ_DIR/$ISO_NAME.collapsed_classification.filtered_lite_classification.txt \
  --iso_gtf $ISOSEQ_SQ_DIR/$ISO_NAME.collapsed_classification.filtered_lite.gtf \
  --ont_gtf $ONT_UNFILTERED_SQ_DIR/$ONT_NAME"_unfiltered_talon_corrected.gtf" \
  --ont_unfil $ONT_UNFILTERED_SQ_DIR/$ONT_NAME"_unfiltered_talon_classification_counts.txt" \
  --ont_fil $ONT_FILTERED_SQ_DIR/$ONT_NAME"_filtered_talon_classification_counts.txt" \
  --cuff $1/1_gffcompare/$NAME.$ONT_NAME"_unfiltered_talon_corrected.gtf.tmap" \
  --o_dir $1/2_common_transcripts \
  --tgenes_ens $2 \
  --sample $3 #&> identify_transcripts.log
    
}

# run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna> <input_dir> 
run_sqanti3(){

    source activate sqanti2_py3
    
    # variable 
    sample=$1
    gtf=$3/$sample.gtf
    #abundance=
    
    # create directory
    mkdir -p $3/3_sqanti3; cd $3/3_sqanti3
    
    # copy STAR output SJ.out.bed files
    SJ_OUT_BED=($(
        for rnaseq in ${RNASEQ_SAMPLES_NAMES[@]}; do
            name=$(find $RNASEQ_MAPPED_DIR -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
            File=$(find $RNASEQ_MAPPED_DIR -name "$name")
            echo $File
        done
        ))
    for file in ${SJ_OUT_BED[@]}; do cp $file $3/3_sqanti3;done
  
    
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


# run_cpat <input_fasta> <output_name> <output_dir>
run_cpat(){
  mkdir -p $3/4_characterise
  mkdir -p $3/4_characterise/CPAT; cd $3/4_characterise/CPAT
  
  source activate sqanti2_py3
  cpat.py -x $HEXAMER -d $LOGITMODEL -g $1 --min-orf=50 --top-orf=50 -o $2 2> $2"_cpat.e"
}


# extract_best_orf <sample> <root_dir>
# extract the best ORF from CPAT for further analysis of ORF predictions for predicted NMD
extract_best_orf(){
  
  # variables
  io_dir=$2/4_characterise/CPAT
  
  cd $io_dir
  extract_fasta_bestorf.py --fa $io_dir/$1".ORF_seqs.fa" --orf $io_dir/$1".ORF_prob.best.tsv" --o_name $1"_bestORF" --o_dir $io_dir &> orfextract.log
}


# convert_gtf_bed12 <input_gtf> 
convert_gtf_bed12(){
  
  # variables 
  output_dir="$(dirname $1)" 
  sample=${1%.gtf} # removes .gtf
  
  source activate sqanti2_py3
  cd $output_dir
  
  gtfToGenePred $1 $sample.genePred
  genePredToBed $sample.genePred > $sample.bed12
  sort -k1,1 -k2,2n $sample.bed12 > $sample"_sorted.bed12"
  rm $sample.genePred $sample.bed12
  # Rscript script.R <input.classfile> <input_output_dir of bed file> <prefix>
  #Rscript $GENERAL/annotate_uscs_tracks.R $SQANTI $WKD $sample
  #bedToBigBed -extraIndex=name $sample $sample"_Modified.bed12" $REFERENCE/mm10.chrom.filtered.sizes $sample.bb
  
  #source deactivate
  #./bedToBigBed -as=bedExample2.as -type=bed9+3 -extraIndex=name $sample"_Modified.bed12" $REFERENCE/mm10.chrom.filtered.sizes $sample.bb
}

# colour_by_abundance <sample> <input_gtf> <abundance_path> <root_dir>
colour_by_abundance(){
  
  mkdir -p $4/4_characterise/bed12Files
  
  convert_gtf_bed12 $2
  bed12=${2%.gtf}_sorted.bed12
  sample="$(basename $2)" 
  outputname=${sample%.gtf}
  echo $outputname
  
  colour_transcripts_by_countandpotential.py \
  --bed $bed12 \
  --cpat $4/4_characterise/CPAT/$1".ORF_prob.best.tsv" \
  --noORF $4/4_characterise/CPAT/$1".no_ORF.txt" \
  --a $3 \
  --o $outputname \
  --dir $4/4_characterise/bed12Files/ \
  --species mouse
  
}


# subset_gene_reference <root_dir>
subset_gene_reference(){
  
  mkdir -p $1/4_characterise/TargetGenesRef
  
  source activate sqanti2_py3
  subset_reference_by_gene.py --r=${GENOME_GTF} --glist ${TGENES[@]} --o $1/4_characterise/TargetGenesRef

}


# run_transdecoder <name> <root_dir>
run_transdecoder(){

  mkdir -p $2/4_characterise/Transdecoder
  cd $2/4_characterise/Transdecoder
  
  source activate nanopore
  TransDecoder.LongOrfs -t $2/4_characterise/CPAT/$1"_bestORF.fasta" &> transdecoder_longorf.log
  
  source activate sqanti2_py3
  hmmsearch --cpu 8 --domtblout pfam.domtblout $PFAM_REF $1"_bestORF.fasta.transdecoder_dir"/longest_orfs.pep &> hmmsearch.log
  
  source activate nanopore
  TransDecoder.Predict -t $2/4_characterise/CPAT/$1"_bestORF.fasta" --retain_pfam_hits pfam.domtblout --no_refine_starts &> transdecoder_predict.log
  sed '/^#/ d' < pfam.domtblout > pfam.domtblout.read
  
}


# full_characterisation <name> <gene> <sq_filter>
full_characterisation(){
  
  mkdir -p $WKD_ROOT/4_characterise/TargetGenes/$3
  mkdir -p $WKD_ROOT/4_characterise/TargetGenes/$3/Log
  
  # variables 
  ref_dir=$WKD_ROOT/4_characterise/TargetGenesRef/
  noISM_path=$WKD_ROOT/3_sqanti3/$1"_ISMrem.classification.txt"
  input_bed=$WKD_ROOT/3_sqanti3/$1.collapsed_corrected_sorted.bed12
  ORF_dir=$WKD_ROOT/4_characterise/CPAT/$1".ORF_prob.best.tsv"  
  output_dir=$WKD_ROOT/4_characterise/TargetGenes/$3/
  
  if [ $3 == "Filtered" ]; then 
    input_gtf=$WKD_ROOT/3_sqanti3/$1.collapsed_classification.filtered_lite.gtf
  elif [ $3 == "Unfiltered" ]; then
    input_gtf=$WKD_ROOT/3_sqanti3/$1.collapsed_corrected.gtf
  else
    echo "Unfiltered/Filtered for 3rd parameter"
  fi
  
  echo "***** $3 dataset"
  echo $input_gtf
  
  source activate sqanti2_py3 
  ficle.py --gene=$2 --ref=$ref_dir \
    --i_bed=$input_bed \
    --i_gtf=$input_gtf \
    --noISM=$noISM_path \
    --orf=$ORF_dir \
    --o_dir=$output_dir &> $output_dir/Log/$2"_characterise.log"
  
}


# run_minimap <output_name> <output_dir> <input_fasta>
run_minimap(){
  
  mkdir -p $2; cd $2
  
  source activate nanopore  
  minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $GENOME_FASTA $3 > $1.sam 2> $1.map.log
  samtools sort -O SAM $1.sam > $1.sorted.sam

  
}


# optimise_collapse <output_name> <input_fasta> <input_map_dir> <max_5_diff> <max_3_diff> 
optimise_collapse(){
  
  output_root_dir="$(dirname $3)"  
  mkdir -p $output_root_dir/$4_5diff_$5_3diff
  cd $output_root_dir/$4_5diff_$5_3diff

  source activate sqanti2_py3
  collapse_isoforms_by_sam.py -c 0.99 -i 0.99 --dun-merge-5-shorter --input $2 -s $3/$1.sorted.sam --max_5_diff=$4 --max_3_diff=5 -o $1 &>> $1.collapse.log

}
