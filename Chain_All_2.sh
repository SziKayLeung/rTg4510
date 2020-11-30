#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:30:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=SQANTI2_test.o
#SBATCH --error=SQANTI2_test.e

source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh

TG_SAMPLES=(K24 L22 M20 O24 P22 Q20 S24 T22     K22     L20     M18     O22     P20     Q18     S22     T20     K20     L18     M24     O20     P18     Q24     S20     T18     K18    L24      M22     O18     P24     Q22     S18     T24)
WT_SAMPLES=(K17 L21 M19 K23 P21 Q19 M21 T21 K21 L19 M17 O21 P19 Q17 S21 T19 K19 L17 M23 O19 P17 Q23 S19 T17 O23 L23 Q21 O17 P23 S23 S17 T23)
SAMPLES_NAMES+=( "${TG_SAMPLES[@]}" "${WT_SAMPLES[@]}" )

################################## All Merged Iso-Seq data

#TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged/TOFU/mm10
KALLISTO=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/Kallisto
#SQANTI2=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/All_Merged
RNASEQ=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/MAPPED/Individual
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019

# run_sqanti2_QC <prefix_sample> <input_tofu_dir> <coverage/genome=mm10_rnqaseq/mm10/hg38_gencode/hg38_chess/ERCC> <input_kallisto_file> <output_dir> <input_rnaseq_dir>
#run_sqanti2_QC All_Merged $TOFU mm10_rnaseq $KALLISTO/All_RNASeq.mod.abundance.tsv $SQANTI2 $RNASEQ

################################## Chained Individual Iso-Seq data
CHAIN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/CHAIN
SQANTI2_POSTCHAIN=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/CHAIN/SQANTI2

#run_sqanti2_QC_postchaining <path_chained.gff> <path_chained_count> <input_kallisto_file> <input_rnaseq_dir> <output_name> <output_dir>
run_sqanti2_QC_postchaining(){
  cd $6

  SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
  CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence

  export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
  python $SQANTI2_DIR/sqanti_qc2.py --version

  # copy STAR output SJ.out.bed files
  SJ_OUT_BED=($(
      for rnaseq in ${SAMPLES_NAMES[@]}; do
          name=$(find $4 -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
          File=$(find $4 -name "$name")
          echo $File
      done
      ))

  for file in ${SJ_OUT_BED[@]}; do cp $file . ;done

  echo "Processing with own RNASeq dataset and genocode.gtf for mouse genome annotation "
  echo "Collapsing with the following bed files"
  ls *SJ.out.bed*


  python $SQANTI2_DIR/sqanti_qc2.py -t 30 \
    --gtf $1 \
    $REFERENCE/gencode.vM22.annotation.gtf  \
    $REFERENCE/mm10.fa \
    --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed \
    --coverage "./*SJ.out.bed" \
    --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt \
    --fl_count $2 \
    --expression $3 \
    --skipORF &> $5.sqanti.qc.log
}

SQANTI3_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/SQANTI3
run_sqanti2_QC_postchaining $CHAIN/all_samples.chained.gff $CHAIN/all_samples.chained_count.txt $KALLISTO/All_RNASeq.mod.abundance.tsv $RNASEQ Individal_PostChain $SQANTI2_POSTCHAIN
cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Tg4510/CHAINED_ANALYIS/mm10/SQANTI2
python $SQANTI2_DIR/sqanti_qc2.py -t 30 ../all_samples.chained.gff $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --expression ../all_samples.chained_count.txt
python $SQANTI3_DIR/sqanti3_qc.py

cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/CHAIN/SQANTI2
python $SQANTI3_DIR/sqanti3_qc.py -t 30 ../all_samples.chained.rep.fq $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --fl_count ../all_samples.chained_count.txt --skipORF
python $SQANTI3_DIR/sqanti3_RulesFilter.py all_samples.chained_classification.txt all_samples.chained_corrected.fasta all_samples.chained_corrected.gtf -a 0.6 -c 3

/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Individual_Samples/CHAIN/SQANTI3
python $SQANTI3_DIR/sqanti3_qc.py -t 30 ../all_samples.chained.rep.fq $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --fl_count ../all_samples.chained_count.txt


INDIVIDUAL_TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Individual/Isoseq/TOFU
L22.collapsed.filtered.abundance.txt L22.collapsed.filtered.gff L22.collapsed.filtered.rep.fq
