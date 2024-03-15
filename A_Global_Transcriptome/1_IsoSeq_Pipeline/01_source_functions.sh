#!/bin/bash
################################################################################################
## Define functions for All_Tg4510_Functions.sh
# 1) run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory>
# 2) run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
# 3) run_REFINE $Sample $Input_LIMA_directory $Output_directory
# 4) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# 5) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
# 6) demux <input path read.stat file> <input path of samples file> <path of output>
# 9) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>

################################################################################################

#************************************* DEFINE VARIABLES
module load Miniconda2/4.3.21
source activate isoseq3

# Listing versions
ccs --version #ccs 5.0.0 (commit v5.0.0)
lima --version #lima 2.0.0 (commit v2.0.0)
isoseq3 --version #isoseq3 3.4.0 (commit v3.4.0)

echo "FASTA SEQUENCE (CLONTECH PRIMERS) FOR NON-MULTIPLEXING"
cat ${FASTA}

################################################################################################
#*************************************  Isoseq3 [Function 1, 2, 3, 4]
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
      time lima $2/$1.ccs.bam ${FASTA} $1.fl.bam --isoseq --dump-clips --dump-removed
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

refine2fasta(){
  input_dir=$1
  mkdir -p $input_dir/fasta
  cd $input_dir/fasta

  Samples=$(echo "${@:2}")
  source activate sqanti2_py3
  for i in ${Samples[@]}; do echo "Processing: $i"; samtools bam2fq $input_dir/$i.flnc.bam | seqtk seq -A > $i.fa; done
  cat *fa* > All_flnc.fasta
  minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 ${GENOME_FASTA} All_flnc.fasta > All.sam 2> All.map.log
  samtools sort -O SAM All.sam > All.sorted.sam
  python ${SEQUENCE}/sam_to_gff3.py All.sorted.sam -i All_flnc.fasta -s ${SPECIES}
}

# run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CLUSTER(){
  source activate isoseq3

  cd $3
  echo "Processing $1 file for cluster"
  if [ -f $1.clustered.bam ]; then
    echo "$1.clustered.bam file already exists; Cluster no need to be processed on Sample $1"
  else
    # cluster <input.flnc.bam> <output.unpolished.bam>
    time isoseq3 cluster $2/$1.flnc.bam $1.clustered.bam --verbose --use-qvs 2> $1.cluster.log
    echo "cluster $1 successful"
    ls $1.clustered*
  fi

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
  isoseq3 --version

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
  source deactivate
}

################################################################################################
#************************************* Post_Isoseq3 (Minimap2, Cupcake, Demultiplex) [Function 5, 6]
# run_pbmm2align <output_name> <clustered_dir> <mapped_dir> 
run_pbmm2align(){

    source activate isoseq3
    echo "Processing Sample $1 for pbmm2 and sort"
    mkdir -p $3; cd $3 #cd to $MAP directory for output
    pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} $2/$1.clustered.hq.fasta $1.bam --log-level TRACE --log-file $1.map.log
    
    source activate nanopore
    samtools view -h $1.bam > $1.sam
    samtools bam2fq $1.bam| seqtk seq -A > $1.fa
    samtools sort -O SAM $1.sam > $1.sorted.sam
}

# filter_alignment <name> <mapped_dir>
filter_alignment(){
  
    source activate nanopore
    
    cd $2
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
}


#  cakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory> 
# Prerequisite: mm10 cage peak
run_map_cupcakecollapse(){

    source activate isoseq3

    samtools --version 
    pbmm2 --version
    isoseq3 collapse --version

    echo "Processing Sample $1 for pbmm2 and sort"
    cd $3 #cd to $MAP directory for output
    pbmm2 align --preset ISOSEQ --sort ${GENOME_FASTA} $2/$1.clustered.hq.fasta $1.bam --log-level TRACE --log-file $1.map.log

    source activate nanopore
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
    mkdir PAF; cd PAF
    source activate nanopore
    htsbox samview -pS $3/$1.sorted.sam > $1.paf
    awk -F'\t' '{if ($6="*") {print $0}}' $1.paf > $1.allread.paf # all reads
    awk -F'\t' '{if ($6=="*") {print $0}}' $1.paf > $1.notread.paf
    awk -F'\t' '{if ($6!="*") {print $0}}' $1.paf > $1.filtered.paf
    awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $1.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $1"_reads_with_alignment_statistics.txt"
    # filter based on alignable length (>0.85) and identity (>0.95)
    awk -F'\t' '{if ($6>=0.85 && $8>=0.95) {print $1}}' $1"_reads_with_alignment_statistics.txt" > $1_filteredreads.txt
    echo "Number of mapped transcripts to genome:"
    wc -l $1.filtered.paf
    echo "Number of ummapped transcripts to genome:"
    wc -l $1.notread.paf
    
    source activate sqanti2
    picard FilterSamReads I=$3/$1.bam O=$3/$1.filtered.bam READ_LIST_FILE=$1_filteredreads.txt FILTER=includeReadList
    
    source activate nanopore
    samtools bam2fq $3/$1.filtered.bam| seqtk seq -A > $3/$1.filtered.fa

    echo "Processing Sample $1 for collapse"
    cd $4 #cd to $TOFU directory for output
    source activate isoseq3
    isoseq3 collapse $3/$1.filtered.bam $1.collapsed.gff --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons --log-level TRACE --log-file $1.collapse.log
  
    echo "Processing Sample $1 for TOFU successful"
}


# demux <name> <refine_dir> <cluster_report> <tofu_dir> 
# run Cupcake_Demultiplex.R, read in read.stat file from cupcake collapse output and count abundance of each sample based on CCS_ID
demux(){
  
    # variables
    refine_dir=$2
    cluster_report=$3
    cluster_dir=$(dirname $3)
    collapse_dir=$4
    
    source activate nanopore
    demux_merged_cluster.py $refine_dir $cluster_report -o $1
    # output = /path/to/cluster_report_dir/name.clustered.demuxed.cluster_report.csv
    
    demux_cupcake_collapse.py $4/$1.collapsed.read_stat.txt ${cluster_dir}/$1.clustered.demuxed.cluster_report.csv --dataset=isoseq -o $1
}



################################################################################################
#************************************* RNASeq & IsoSeq [Function 9]
# mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
# Aim: Merge forward and reverse of all RNASeq filtered fastq, for further downstream alignment to IsoSeq Tofu output using Kallisto
mouse_merge_fastq(){
    R1_READS=($(
        for i in ${RNASEQ_SAMPLES_NAMES[@]}; do
            F_name=$(find $1 -name "*fastq.filtered" -exec basename \{} \; | grep ^$i | grep "R1" )
            F_File=$(find $1 -name "$F_name")
            echo $F_File
        done
    ))

    R2_READS=($(
        for i in ${RNASEQ_SAMPLES_NAMES[@]}; do
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
#************************************* SQANTI2 [Function 10]
# run_sqanti3 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/rnaseq/lncrna>
run_sqanti3(){

    source activate sqanti2_py3

    # copy STAR output SJ.out.bed files
    SJ_OUT_BED=($(
        for rnaseq in ${RNASEQ_SQ_INPUT[@]}; do
            name=$(find $4 -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
            File=$(find $4 -name "$name")
            echo $File
        done
        ))
    for file in ${SJ_OUT_BED[@]}; do cp $file $7 ;done

    # create tab separated file from kallisto
    # Rscript script.R <input.file_fromkallisto> <output.file_intosqanti>
    Rscript ${MODKAL} $5 $7/$1".mod2.abundance.tsv"
    kallisto_expfile=$7/$1".mod2.abundance.tsv"
    echo "Using $kallisto_expfile"

    # prepare sqanti
    mkdir -p $7; cd $7
    export PYTHONPATH=$PYTHONPATH:${SEQUENCE}
    python ${SQANTI3_DIR}/sqanti3_qc.py -v

    # sqanti qc
    echo "Processing Sample $1 for SQANTI2 QC"
    
    # no kalliso file
    if [ $8 == "rnaseq" ]; then
      python ${SQANTI3_DIR}/sqanti3_qc.py -t 30 $3/$2 ${GENOME_GTF} ${GENOME_FASTA} --CAGE_peak ${CAGE_PEAK} --coverage "./*SJ.out.bed" --polyA_motif_list ${POLYA} --genename --isoAnnotLite --gff3 ${GFF3} --report skip &> $1.sqanti.qc.log

    elif [ $8 == "lncrna" ]; then
      echo "Processing with lncRNA.gtf for genome annotation "
      python ${SQANTI3_DIR}/sqanti3_qc.py -t 30 $3/$2 ${GENOME_LNCRNA_GTF} ${GENOME_FASTA} --CAGE_peak ${CAGE_PEAK} --coverage "./*SJ.out.bed" --polyA_motif_list ${POLYA} --skipORF --fl_count $6 --genename --isoAnnotLite --gff3 ${GFF3} --report skip &> $1.sqanti.qc.log

    elif [ $8 == "genome" ]; then
      echo "Processing with gencode.vM22.annotation.gtf for genome annotation "
      python ${SQANTI3_DIR}/sqanti3_qc.py -t 30 $3/$2 ${GENOME_GTF} ${GENOME_FASTA} --CAGE_peak ${CAGE_PEAK} -c ./"*SJ.out.bed" --polyA_motif_list ${POLYA} --expression $kallisto_expfile --fl_count $6 --genename --isoAnnotLite --gff3 ${GFF3} --report skip &> $1.sqanti.qc.log
    
    elif [ $8 == "basic" ]; then
      echo "Processing with gencode.vM22.annotation.gtf for genome annotation and basic"
      python ${SQANTI3_DIR}/sqanti3_qc.py -t 30 $3/$2 ${GENOME_GTF} ${GENOME_FASTA} --CAGE_peak ${CAGE_PEAK} --polyA_motif_list ${POLYA} --genename --isoAnnotLite --gff3 ${GFF3} --report skip &> $1.sqanti.qc.log
      
    elif [ $8 == "nornaseq" ]; then
      echo "Processing with gencode.vM22.annotation.gtf for genome annotation and no RNASeq Junction"
      python ${SQANTI3_DIR}/sqanti3_qc.py -t 30 $3/$2 ${GENOME_GTF} ${GENOME_FASTA} --CAGE_peak ${CAGE_PEAK} --polyA_motif_list ${POLYA} --fl_count $6 --genename --isoAnnotLite --gff3 ${GFF3} --report skip &> $1.sqanti.qc.log


    else
      echo "8th argument required"
    fi

    echo "Processing Sample $1 for SQANTI2 filter"
    filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
    python $SQANTI3_DIR/sqanti3_filter.py rules $1"_classification.txt" --faa=$1"_corrected.fasta" --gtf=$1"_corrected.gtf" -j=${filteringJson} --skip_report &> $1.sqanti.filter.log

    # remove temp SJ.out bed files
    rm *SJ.out.bed
    source deactivate
}

