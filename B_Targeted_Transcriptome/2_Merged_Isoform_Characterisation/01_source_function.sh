#run_gff_compare <output_name> <ref_gtf> <2nd_gtf> <root_dir>
run_gff_compare(){
  
  mkdir -p $WKD_ROOT/1_gffcompare; cd $WKD_ROOT/1_gffcompare
  cp $2 . 
  cp $3 . 
  
  source activate sqanti2_py3 
  gffcompare -r $2 $3 -o $1
  
}

# identify_common_transcripts <root_dir>
# Custom script to identify_common_targeted_transcripts for merged output
identify_common_transcripts(){
  
  mkdir -p $1/2_common_transcripts; cd $1/2_common_transcripts
  
  source activate sqanti2_py3
  python $COMMONSC \
  --iso $ISOSEQ_SQ_DIR/$ISO_NAME.collapsed_classification.filtered_lite_classification.txt \
  --iso_gtf $ISOSEQ_SQ_DIR/$ISO_NAME.collapsed_classification.filtered_lite.gtf \
  --ont_gtf $ONT_UNFILTERED_SQ_DIR/$ONT_NAME"_unfiltered_talon_corrected.gtf" \
  --ont_1 $ONT_UNFILTERED_SQ_DIR/$ONT_NAME"_unfiltered_talon_classification.txt" \
  --ont_2 $ONT_FILTERED_SQ_DIR/$ONT_NAME"_filtered_talon_classification.txt" \
  --a_ont $ONT_UNFILTERED_DIR/$ONT_NAME"_unfiltered_talon_abundance.tsv" \
  --cuff $1/1_gffcompare/$NAME.$ONT_NAME"_unfiltered_talon_corrected.gtf.tmap" \
  --o_dir 2_common_transcripts &> identify_transcripts.log
    
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
  mkdir -p $3/4_characterise; cd $3/4_characterise/CPAT
  
  source activate sqanti2_py3
  cpat.py -x $HEXAMER -d $LOGITMODEL -g $1 --min-orf=50 --top-orf=50 -o $2 2> $2"_cpat.e"
}


# extract_best_orf <sample> <root_dir>
# extract the best ORF from CPAT for further analysis of ORF predictions for predicted NMD
extract_best_orf(){
  
  # variables
  io_dir=$2/4_characterise/CPAT
  
  cd $io_dir
  python $BESTORF --fa $io_dir/$1".ORF_seqs.fa" --orf $io_dir/$1".ORF_prob.best.tsv" --o_name $1"_bestORF" --o_dir $io_dir &> orfextract.log
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

# colour_by_abundance <sample> <input_gtf> <root_dir>
colour_by_abundance(){
  
  convert_gtf_bed12 $2
  bed12=${2%.gtf}_sorted.bed12
  sample="$(basename $2)" 
  output_bed12=${sample%.gtf}_sorted_coloured.bed12
  
  python $COLOURTRANS \
  --bed $bed12 \
  --cpat $3/4_characterise/CPAT/$1".ORF_prob.best.tsv" \
  --noORF $3/4_characterise/CPAT/$1".no_ORF.txt" \
  --a $3/2_common_transcripts/Final_Merged_Abundance.csv \
  --o $3/4_characterise/$output_bed12
  
}


# subset_gene_reference 
subset_gene_reference(){
  
  mkdir -p $WKD_ROOT/4_characterise/TargetGenesRef
  
  source activate sqanti2_py3
  python $TGENEPREP --r=$GENOME_GTF --glist ${TGENES[@]} --o $WKD_ROOT/4_characterise/TargetGenesRef

}


# run_transdecoder <name> <root_dir>
run_transdecoder(){
  source deactivate
  
  mkdir -p $2/4_characterise/Transdecoder; cd $2/4_characterise/Transdecoder
  
  TransDecoder.LongOrfs -t $2/4_characterise/CPAT/$1"_bestORF.fasta" &> transdecoder_longorf.log
  
  source activate sqanti2_py3
  hmmscan --cpu 8 --domtblout pfam.domtblout $PFAM_REF $1"_bestORF.fasta.transdecoder_dir"/longest_orfs.pep &> hmmscan.log
  
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
  python $FICLE --gene=$2 --ref=$ref_dir \
  --i_bed=$input_bed \
  --i_gtf=$input_gtf \
  --noISM=$noISM_path \
  --orf=$ORF_dir \
  --o_dir=$output_dir &> $output_dir/Log/$2"_characterise.log"
  
}