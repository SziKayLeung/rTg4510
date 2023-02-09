# run_kallisto_1sample <input_RNASEQ_rawdir> <sample> <input_ref_name_idx> <output_dir>
run_kallisto_1sample(){
    source activate sqanti2_py3

    # variables
    input_dir=$1
    sample=$2
    input_ref_name=$3
    output_dir=$4
    
    # cd to $J20/Tg4510_input_directory
    cd $input_dir
    echo "Finding reads for Sample $sample for STAR"
    # find reverse and forward file, trimming "./" therefore only printing file name
    F=$(find . | grep "fastq" | grep $sample | grep "R1" | sed 's|^./||' )
    R=$(find . | grep "fastq" | grep $sample | grep "R2" | sed 's|^./||' )
    # save path directory of files as variable for later mapping
	  R1_READS=$(realpath $F)
    R2_READS=$(realpath $R)

    echo "Processing Kallisto for $sample"
    echo $R1_READS
    echo $R2_READS

    cd $output_dir
    #kallisto version
    time kallisto quant -i $input_ref_name --rf-stranded $R1_READS $R2_READS -o $sample 2> $2"_Kallisto.quant.log"
 
    source deactivate
}

# generate_rnaseq_counts <input_dir>
# input_dir = directory containing kallisto output files from alignment of RNA-Seq to Iso-Seq annotation file
generate_rnaseq_counts(){
  
  source activate nanopore
  # Rscript script.R <input.dir> <output.file> <type=Whole/Targeted/WholeTargeted> <targeted.class.files>
  Rscript $RNASEQCOUNT $1 $NAME"_RNASeq.expression.txt" Whole NA
  
}