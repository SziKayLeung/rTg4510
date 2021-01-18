#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Individual_Part5.o
#SBATCH --error=Individual_Part5.e

# 04/01/2021: post tofu on all samples, demultiplex using readstats and SQANTI

# File directories
SQANTI=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/SQANTIDEMUX
SQANTIAll=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/SQANTIAll
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Samples/ERCC/TOFU

# run_sqanti2 <sample_name> <working_directory> <collapsed.gtf> <count.txt> <mm10/ERCC>
run_sqanti2(){
    module load Miniconda2
    source activate sqanti2_py3
    SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    CAGE_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CAGE

    echo "Processing Sample $1 for SQANTI2 QC"
    cd $2

    export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
    python $SQANTI2_DIR/sqanti_qc2.py -v

    if [ $5 == "mm10" ]; then
      if python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
          --gtf $3 \
          $REFERENCE/gencode.vM22.annotation.gtf  \
          $REFERENCE/mm10.fa \
          --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed \
          --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt \
          --fl_count $4 \
          2> $1.sqanti.qc.log; then
          echo "Processing Sample $1 for SQANTI2 successful"
      else
          echo "Processing Sample $1 for SQANTI2 failed"
      fi
     elif [ $5 == "ERCC" ]; then
       if python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
          --gtf $3 \
          $REFERENCE/ERCC/ERCC92.gtf  \
          $REFERENCE/ERCC/ERCC92.fa \
          --fl_count $4 \
          2> $1.sqanti.qc.log; then
          echo "Processing Sample $1 for SQANTI2 successful"
      else
          echo "Processing Sample $1 for SQANTI2 failed"
      fi
    else
        echo "ERROR: Require 5th argument to be either: mm10 or ERCC"
        return
    fi

    mkdir Unfiltered_PNG; mv *png* Unfiltered_PNG/

    if python $SQANTI2_DIR/sqanti_filter2.py \
        $1"_classification.txt" \
        $1"_corrected.fasta"  \
        $1"_corrected.gtf" \
        -a 0.6 -c 3 2> $1.sqanti.filter.log; then
        echo "Processing Sample $1 for SQANTI2 filter successful"
    else
        echo "Processing Sample $1 for SQANTI2 filter failed"
    fi

    mkdir filtered_PNG; mv *png* filtered_PNG/
    source deactivate

}

run_sqanti2 All_Merged.collapsed.filtered $SQANTI $TOFU/All_Merged.collapsed.filtered.gff $SQANTI/Mouse.Demultiplexed_Abundance.txt ERCC
run_sqanti2 All_Merged.collapsed.filtered $SQANTIAll $TOFU/All_Merged.collapsed.filtered.gff $TOFU/All_Merged.collapsed.filtered.abundance.txt ERCC
