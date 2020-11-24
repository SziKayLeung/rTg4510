#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=10:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# 17/10/2019: run SQANTI2 using count.chain txt files as input 

#************************************* DEFINE GLOBAL VARIABLES
module load Miniconda2
source activate sqanti2

SQANTI2_output_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/CHAIN/SQANTI2
CHAIN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/CHAIN
STAR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/MAPPED/Individual
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/TOFU
SQANTI2_filter_output_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq3.2.1/CHAIN/SQANTI2

#************************************* SQANTI2 FUNCTION
run_sqanti2_QC(){

    SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    CAGE_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CAGE

    cd $SQANTI2_output_dir
        
    for i in "$@"; do 
        echo "Processing Sample $i for SQANTI2"
        export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
        python $SQANTI2_DIR/sqanti_qc2.py -v

        echo "Finding RNA-Seq data for Sample $i for input to SQANTI2"
        Samples_SJ_out_bed=$(find $STAR -name "*SJ.out.bed" | grep $i )
        echo ${Samples_SJ_out_bed[@]}

        if python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
        $TOFU/$i.collapsed.filtered.rep.fa \
        $REFERENCE/gencode.vM22.annotation.gtf  \
        $REFERENCE/mm10.fa \
        --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed \
        --coverage "$Samples_SJ_out_bed" \
        --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt \
        --fl_count $CHAIN/all_samples.chained_count.txt \
        &> $i.sqanti.qc.log; then
            echo "Processing Sample $i for SQANTI2 successful" # convert SQANTI2 output sam file to output bam file for Alignqc input 
            samtools view -S -b $i.collapsed.filtered.rep.renamed_corrected.sam > $i.collapsed.filtered.rep.renamed_corrected.bam 
        else 
            echo "Processing Sample $i for SQANTI2 failed" 
        fi   
    done
}

run_sqanti2_Filter(){

    cd $SQANTI2_filter_output_dir  
  
    for i in "$@"; do 
    echo "Processing Sample $i for SQANTI2 filter"
    export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
    if python $SQANTI2_DIR/sqanti_filter2.py \
        $SQANTI2_output_dir/$i.collapsed.filtered.rep_classification.txt \
        $SQANTI2_output_dir/$i.collapsed.filtered.rep.renamed_corrected.fasta \
        $SQANTI2_output_dir/$i.collapsed.filtered.rep.renamed_corrected.sam \
        $SQANTI2_output_dir/$i.collapsed.filtered.rep.renamed_corrected.gtf \
        -a 0.6 -c 3; then
            echo "Processing Sample $i for SQANTI2 filter successful"
        else
            echo "Processing Sample $i for SQANTI2 filter failed"
        fi
    done
}

run_sqanti2_QC Q21
run_sqanti2_Filter Q21