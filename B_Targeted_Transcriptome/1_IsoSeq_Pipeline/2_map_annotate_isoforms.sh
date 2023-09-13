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
#SBATCH --output=2_map_annotate_isoforms.o
#SBATCH --error=2_map_annotate_isoforms.e

# 11/05/2021: Run pipeline from refine to sqanti and tama (output: Targeted_ADBDR_Part2.o)
# 03/11/2021: Run kallisto with ADBDR RNA-Seq (30 samples), SQANTI3 with RNA-Seq support and TAMA filter (output: Targeted_ADBDR_Part2b.o)

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_isoseq_processing

export samplename=AllMouseTargeted
mkdir -p $WKD_ROOT/6_minimap $WKD_ROOT/7a_targeted_only

##-------------------------------------------------------------------------
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch (1_run_isoseq3.sh) for all batched runs"

## 4) run_targeted_REFINE 
run_targeted_REFINE 

# 6) run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
for i in ${ALL_SAMPLES_NAMES[@]}; do run_CLUSTER $i; done

##-------------------------------------------------------------------------
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
# merging_at_refine <output_name> <output_root_dir> <samples.....>
merging_at_refine $NAME $WKD_ROOT $WKD_ROOT ${ALL_SAMPLES_NAMES[@]}

# run_pbmm2align <output_name> <clustered_dir> <mapped_dir> 
for i in ${ALL_SAMPLES_NAMES[@]}; do run_pbmm2align $i $WKD_ROOT/4_cluster $WKD_ROOT/6_minimap; done

# run_target_rate

# filter_alignment <name> <mapped_dir>
for i in ${ALL_SAMPLES_NAMES[@]}; do filter_alignment $i $WKD_ROOT/6_minimap; done

# merge_run_isoseqcollapse <name> <mapped_dir> <outputdir> <samples....>  
merge_run_isoseqcollapse ${samplename} $WKD_ROOT/6_minimap $WKD_ROOT/7a_targeted_only ${ALL_SAMPLES_NAMES[@]}

# demux <name> <refine_dir> <cluster_report> <tofu_dir> 
demux ${samplename} $WKD_ROOT/3_refine $WKD_ROOT/5_merged_cluster/AllMouseTargeted.clustered.cluster_report.csv $WKD_ROOT/7a_targeted_only

# sqanti3
echo "Running SQANTI3..."
mkdir -p $WKD_ROOT/7a_targeted_only/8_sqanti3; cd $WKD_ROOT/7a_targeted_only/8_sqanti3
source activate sqanti2_py3
python $SQANTI3_DIR/sqanti3_qc.py $WKD_ROOT/7a_targeted_only/${samplename}.gff \
  $GENOME_GTF $GENOME_FASTA -t 30 --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report skip --fl_count $WKD_ROOT/7a_targeted_only/${samplename}_fl_count.csv &> ${samplename}_sqanti_qc.log

# sqanti3 filter 
filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_classification.txt \
  --faa=${samplename}_corrected.fasta \
  --gtf=${samplename}_corrected.gtf \
  -j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log



##-------------------------------------------------------------------------
# find transgene

find_transgene(){
  # create a file listing all the readID in cluster report (output = pre_cluster_read.csv)
  mkdir -p $WKD_ROOT/10_characterise/transgene
  cd ${WKD_ROOT}/4_cluster/
  for i in *cluster_report.csv*; do wc -l $i; done > pre_cluster_read.csv
  sed -i 's/ \+/,/g' pre_cluster_read.csv 
  sed -i '1 i num_reads,file' pre_cluster_read.csv
  mv pre_cluster_read.csv $WKD_ROOT/10_characterise/transgene
  
  
  # grep transgene sequences in clustered .fasta
  name=(hmapt1 hmapt2 mmapt1)
  seq=($hMAPT_1 $hMAPT_2 $mMAPT_1)
  for i in {0..2}; do
    echo ${name[$i]}
    echo ${seq[$i]}
    search_fasta_by_sequence.py --fasta=${WKD_ROOT}/4_cluster --i=hq.fa \
      --seq=${seq[$i]} -o=${name[$i]} \
      -d=$WKD_ROOT/10_characterise/transgene
  done
}