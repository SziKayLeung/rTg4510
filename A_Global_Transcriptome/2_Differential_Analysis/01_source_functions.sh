#!/bin/bash


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

      # Parameters according to https://github.com/jennylsmith/Isoseq3_workflow/blob/master/shell_scripts/8__STAR_Junctions_ShortReads.sh
      # 2-pass mode alignment with chimeric read detection
      # at least 25 bp of one read of a pair maps to a different loci than the rest of the read pair
      # require 20 pb on either side of chimeric junction
      # include chimeric reads in the output BAM
      # don't include chimeras that map to reference genome contig with Ns
      # --outSAMtype BAM SortedByCoordinate, output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.
      # note STAR indexed with genocde vM22 primary assembly annotation gtf therefore no need to include into command line, otherwise create additional folders
        STAR --runMode alignReads --runThreadN 32 --genomeDir ${STAR_REF_DIR} \
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
        for i in ${RNASEQ_SQ_INPUT[@]}; do
            F_name=$(find $1 -name "*fastq.filtered" -exec basename \{} \; | grep ^$i | grep "R1" )
            F_File=$(find $1 -name "$F_name")
            echo $F_File
        done
    ))

    R2_READS=($(
        for i in ${RNASEQ_SQ_INPUT[@]}; do
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

# run_kallisto_1sample <input_RNASEQ_rawdir> <sample> <input_ref_name_idx> <output_dir>
run_kallisto_1sample(){
    source activate sqanti2

    # variables
    input_dir=$1
    sample=$2
    input_ref_name=$3
    output_dir=$4


    R1_READS=($(
          F_name=$(find $input_dir -name "*fastq.filtered" -exec basename \{} \; | grep ^$sample | grep "R1" )
          F_File=$(find $input_dir -name "$F_name")
          echo $F_File
    ))


    R2_READS=($(
          R_name=$(find $input_dir -name "*fastq.filtered" -exec basename \{} \; | grep ^$sample | grep "R2" )
          R_File=$(find $input_dir -name "$R_name")
          echo $R_File
    ))


    echo "Processing Kallisto for $sample"
    echo $R1_READS
    echo $R2_READS

    cd $output_dir
    #kallisto version
    time kallisto quant -i $input_ref_name --rf-stranded $R1_READS $R2_READS -o $sample 2> $sample"_Kallisto.quant.log"

    source deactivate
}

################################################################################################
#************************************* Prepare Input files for tappAS [Function 7]

# counts_subset_4tappas <input_class> <output_class> <type_genes>
# regenerate expression matrix for tappas from long reads FL read counts
counts_subset_4tappas(){
  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

  # variable
  input_class=$1
  output_class=$2
  type_genes=$3

  # Rscript script.R <input.classfile> <output.classfile> <type>
  # type = AD or Non-AD
  source activate sqanti2_py3
  Rscript $GENERALFUNC/Counts_Subset.R $input_class $output_class $type_genes
}

################################################################################################
#************************************* Alternative Splicing [Function 8]

# run_suppa2 <input_gtf> <input_class> <output_dir> <output_name>
# input_gtf: sqanti tama filtered gtf
# input_class: sqanti tama filtered classification file
run_suppa2(){
  # variable
  input_gtf=$1
  input_class=$2
  output_dir=$3
  output_name=$4

  echo "#******************* $output_name"
  echo "Processing gtf: $input_gtf"
  echo "Processing classification file: $input_class"

  FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

  source activate sqanti2_py3
  cd $output_dir
  suppa.py generateEvents -i $input_gtf -o $output_name --pool-genes -f ioe -e SE MX FL SS RI &>> $output_name"_suppa2.log"
  python $FUNCTIONS/Suppa2_output_mod_updated.py $output_name $output_dir $input_class &>> $output_name"_suppa2_mod.log"
}

################################################################################################
#************************************* Whole + Targeted Transcriptome
# convertgtf2bed <name> <input_gtf> <working_directory>
convertgtf2bed(){
  # variable
  name=$1
  input_gtf=$2
  wkd=$3

  cd $wkd

  source activate sqanti2
  # convert gtf to bed
  gtfToGenePred $input_gtf $name.genepred
  genePredToBed $name.genepred $name.bed12

  # change bed into correct format for tama
  awk -F'\t' '{print $1,$2,$3,$4,"40",$6,$7,$8,"255,0,0",$10,$11,$12}' $name.bed12 | sed s/" "/"\t"/g|sed s/",\t"/"\t"/g|sed s/",$"/""/g > $name.mod.bed12
  TAMA_MERGE_prepare=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/TAMA/TAMA_Merge_Prepare.R
  # Rscript script.R <name>.bed12 <input_dir>
  Rscript $TAMA_MERGE_prepare $name".mod" $wkd
}

# run_tamamerge <wkd> <wholebed> <targetedbed> <outputname>
run_tamamerge(){
  # variable
  wkd=$1
  wholebed=$2
  targetedbed=$3
  outputname=$4

  # prepare file list
  echo "Merging the following files for TAMA"
  echo "$wholebed:no_cap:1,1,1:Whole"|sed s/":"/"\t"/g>file.list
  echo "$targetedbed:no_cap:1,1,1:Targeted"|sed s/":"/"\t"/g>>file.list
  cat file.list

  # run tama merge using the file list
  cd $wkd
  source activate sqanti2
  tamasoftware=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama
  python $tamasoftware/tama_merge.py -f file.list -p $outputname &> tama_output.log

  # variable
  wkd=$1
  wholebed=$2
  targetedbed=$3
  outputname=$4
  # convert from bed to gtf
  # bedToGenePred $outputname.bed $outputname.genepred
  # genePredToGtf $outputname.genepred $outputname.gtf
  source activate sqanti2
  #bedToGenePred $outputname.bed /dev/stdout | genePredToGtf file /dev/stdin $outputname.gtf
  tamasoftware=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama/tama_go/format_converter/
  python $tamasoftware/tama_convert_bed_gtf_ensembl_no_cds.py $outputname.bed $outputname.gtf
}
