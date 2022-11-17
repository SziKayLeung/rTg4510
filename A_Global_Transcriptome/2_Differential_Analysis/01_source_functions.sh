#!/bin/bash
################################################################################################
## Define functions for Diff_Whole.sh
## SQANTI3 & TAMA Filter
# 1) run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/rnaseq/lncrna>
# 2) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# 3) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <sqanti_output_junc> <output_prefix_name> <output_dir>
## RNA-Seq Expression Matrix on Iso-Seq scaffold
# 4) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
# 5) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# 6) run_kallisto_1sample <input_RNASEQ_rawdir> <sample> <input_ref_name_idx> <output_dir>
## Prepare Input files for tappAS
# 7) counts_subset_4tappas <input_class> <output_class> <type_genes>


################################################################################################
#************************************* SQANTI3 & TAMA Filter [Function 1 - 3]
# run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/rnaseq/lncrna>
run_sqanti3(){

    source activate sqanti2_py3

    # variables
    SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    SQANTI3_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/SQANTI3
    CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    CAGE_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CAGE
    TAPPAS_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/TAPPAS
    GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

    # copy STAR output SJ.out.bed files
    SJ_OUT_BED=($(
        for rnaseq in ${SAMPLES_NAMES[@]}; do
            name=$(find $4 -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
            File=$(find $4 -name "$name")
            echo $File
        done
        ))
    for file in ${SJ_OUT_BED[@]}; do cp $file $7 ;done

    # create tab separated file from kallisto
    # Rscript script.R <input.file_fromkallisto> <output.file_intosqanti>
    Rscript $GENERALFUNC/TabSeparated_Kallisto.R $5 $7/$1".mod2.abundance.tsv"
    kallisto_expfile=$7/$1".mod2.abundance.tsv"
    echo "Using $kallisto_expfile"

    # prepare sqanti
    cd $7
    export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
    python $SQANTI3_DIR/sqanti3_qc.py -v

    # sqanti qc
    echo "Processing Sample $1 for SQANTI2 QC"

    # no kalliso file
    if [ $8 == "rnaseq" ]; then
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --cage_peak $SQANTI3_DIR/data/ref_TSS_annotation/mouse.refTSS_v3.1.mm10.bed --coverage "./*SJ.out.bed" --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --genename --isoAnnotLite --gff3 $TAPPAS_dir/Mus_musculus_GRCm38_Ensembl_86.gff3 --report pdf &> $1.sqanti.qc.log

    elif [ $8 == "lncrna" ]; then
      echo "Processing with lncRNA.gtf for genome annotation "
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM25.long_noncoding_RNAs.gtf $REFERENCE/mm10.fa --cage_peak $SQANTI3_DIR/data/ref_TSS_annotation/mouse.refTSS_v3.1.mm10.bed --coverage "./*SJ.out.bed" --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --skipORF --fl_count $6 --genename --isoAnnotLite --gff3 $TAPPAS_dir/Mus_musculus_GRCm38_Ensembl_86.gff3 --report pdf &> $1.sqanti.qc.log

    elif [ $8 == "genome" ]; then
      echo "Processing with gencode.vM22.annotation.gtf for genome annotation "
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --cage_peak $SQANTI3_DIR/data/ref_TSS_annotation/mouse.refTSS_v3.1.mm10.bed -c ./"*SJ.out.bed" --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --expression $kallisto_expfile --fl_count $6 --genename --isoAnnotLite --gff3 $TAPPAS_dir/Mus_musculus_GRCm38_Ensembl_86.gff3 --report pdf &> $1.sqanti.qc.log
    
    elif [ $8 == "basic" ]; then
      echo "Processing with gencode.vM22.annotation.gtf for genome annotation and basic"
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --cage_peak $SQANTI3_DIR/data/ref_TSS_annotation/mouse.refTSS_v3.1.mm10.bed --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --genename --isoAnnotLite --gff3 $TAPPAS_dir/Mus_musculus_GRCm38_Ensembl_86.gff3 --report pdf &> $1.sqanti.qc.log
      
    elif [ $8 == "nornaseq" ]; then
      echo "Processing with gencode.vM22.annotation.gtf for genome annotation and no RNASeq Junction"
      python $SQANTI3_DIR/sqanti3_qc.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --cage_peak $SQANTI3_DIR/data/ref_TSS_annotation/mouse.refTSS_v3.1.mm10.bed --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --fl_count $6 --genename --isoAnnotLite --gff3 $TAPPAS_dir/Mus_musculus_GRCm38_Ensembl_86.gff3 --report pdf &> $1.sqanti.qc.log


    else
      echo "8th argument required"
    fi

    echo "Processing Sample $1 for SQANTI2 filter"
    python $SQANTI3_DIR/sqanti3_RulesFilter.py $1"_classification.txt" $1"_corrected.fasta" $1"_corrected.gtf" -a 0.6 -c 3 &> $1.sqanti.filter.log

    # remove temp SJ.out bed files
    rm *SJ.out.bed
    source deactivate
}


# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# remove short fragments from post tofu
# Prerequisite: Require TAMA_prepare.R to change column 4 of bed12 to gene_name: transcript_name for correct file TAMA format
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

  # reinsert the quotation marks around the transcript id, gene id etc as required for recognition by SQANTI2
  # note, these are removed after processing in R
  sed 's/transcript_id \([^;]\+\)/transcript_id \"\1\"/g' $7"_sqantitamafiltered.classification.gtf" | sed 's/gene_id \([^;]\+\)/gene_id \"\1\"/g' | sed 's/gene_name \([^;]\+\)/gene_name \"\1\"/g' | sed 's/ref_gene_id \([^;]\+\)/ref_gene_id \"\1\"/g' > $7"_sqantitamafiltered.final.classification.gtf"
}

collapse_filter(){
  source activate sqanti2_py3
  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
# Rscript .R <path/tama_filtered_output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_junc_txt> <output_prefix_name> <output_dir>
  Rscript $GENERALFUNC/sqanti_classgtfsubset.R $1 $2 $3 $4 $6 $7 $8 &> $8/$7"_subset.classgtf.log"
  python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $2/$5 $1 $8/$7"_sqantisubset_classification.fasta" &> $8/$7"_subset.fasta.log"

  # reinsert the quotation marks around the transcript id, gene id etc as required for recognition by SQANTI2
  # note, these are removed after processing in R
  cd $8
  sed 's/transcript_id \([^;]\+\)/transcript_id \"\1\"/g' $7"_sqantisubset.classification.gtf" | sed 's/gene_id \([^;]\+\)/gene_id \"\1\"/g' | sed 's/gene_name \([^;]\+\)/gene_name \"\1\"/g' | sed 's/ref_gene_id \([^;]\+\)/ref_gene_id \"\1\"/g' > $7"_sqantisubset.final.classification.gtf"
}



################################################################################################
#************************************* RNA-Seq Expression Matrix on Iso-Seq scaffold [Function 4, 5, 6]
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


# run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# Aim: Align merged RNASeq fastq files to IsoSeq Tofu output fasta using Kallisto
# Prerequisite:
    # run mouse_merge_fastq/any equivalent merging of all RNASeq fastq files (note sample_prefix_output name is the same)
# Output: <output_prefix>.mod.abundance.tsv for input into Sqanti_qc.py (not modified script to take in space delimited rather than tab file)
run_kallisto(){
    #source activate nanopore
    #check_strandedness --gtf $REFERENCE/gencode.vM22.annotation.gtf --transcripts $REFERENCE/mm10transcriptome.fa --reads_1 $RNASeq_Filtered/Tg4510_filtered/K17/K17_S1_R1_001.fastq.filtered --reads_2 $RNASeq_Filtered/Tg4510_filtered/K17/K17_S1_R2_001.fastq.filtered

    source activate sqanti2

    echo "Processing Kallisto for $1"
    cd $4
    kallisto version
    time kallisto index -i $1_Kallisto.idx $2 2> $1_Kallisto.index.log
    time kallisto quant -i $4/$1_Kallisto.idx --rf-stranded $3/$1_R1.fq $3/$1_R2.fq -o $4 2> $1_Kallisto.quant.log
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
