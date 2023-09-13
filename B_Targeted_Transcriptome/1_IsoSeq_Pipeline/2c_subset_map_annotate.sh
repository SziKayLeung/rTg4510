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
#SBATCH --output=2c_subset_map_annotate.o
#SBATCH --error=2c_subset_map_annotate.3

# 07/06/2023: subset matched samples from whole and targeted Iso-Seq
  # rerun cluster with only matched samples
  # align using pbmm2
  # isoseq3 collapse
  # sqanti3 QC and filter
  
##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/rTg4510_isoseq.config
source $SC_ROOT/B_Targeted_Transcriptome/1_IsoSeq_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/assist_isoseq_processing
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing

export samplename=MatchedMouse
mkdir -p $WKD_ROOT/7b_matched_only

export GLOBAL_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/A_IsoSeq_Whole

##-------------------------------------------------------------------------
run_cupcake(){
  # merge flnc files of matching samples from whole and targeted Iso-Seq dataset
  echo "Merging flnc bam files..."
  globIsoseq_flnc=$(for i in ${SUBSET_SAMPLES_NAMES[@]}; do echo $(find ${GLOBAL_ROOT}/1_isoseq3/3_refine -name "*.flnc.bam" | grep "$i" ); done)
  targIsoseq_flnc=$(for i in ${SUBSET_SAMPLES_NAMES[@]}; do echo $(find $WKD_ROOT/3_refine -name "*.flnc.bam" | grep "$i" ); done)
  all_flnc_bams=(${globIsoseq_flnc[@]} ${targIsoseq_flnc[@]})
  mkdir -p $WKD_ROOT/7b_matched_only/8a_merged_cluster; cd $WKD_ROOT/7b_matched_only/8a_merged_cluster
  printf '%s\n' "${all_flnc_bams[@]}" > $samplename.flnc.fofn
  cat $samplename.flnc.fofn
  
  # cluster only those samples
  source activate isoseq3
  isoseq3 cluster $samplename.flnc.fofn $samplename.clustered.bam --verbose --use-qvs
  gunzip $samplename.clustered.hq.fasta.gz
  
  # alignment of clustered samples 
  mkdir -p $WKD_ROOT/7b_matched_only/8b_pbmm2
  run_pbmm2align $samplename $WKD_ROOT/7b_matched_only/8a_merged_cluster $WKD_ROOT/7b_matched_only/8b_pbmm2
  filter_alignment $samplename $WKD_ROOT/7b_matched_only/8b_pbmm2
  
  # collapse
  source activate isoseq3
  mkdir -p $WKD_ROOT/7b_matched_only/8c_collapse; cd $WKD_ROOT/7b_matched_only/8c_collapse
  isoseq3 collapse $WKD_ROOT/7b_matched_only/8b_pbmm2/${samplename}.filtered.bam ${samplename}.gff --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
      --log-level TRACE --log-file ${samplename}.log
      
  # demux 
  globIsoseq_cluster=$(for i in ${SUBSET_SAMPLES_NAMES[@]}; do echo $(find ${GLOBAL_ROOT}/1_isoseq3/3_refine -name "*flnc.report.csv" | grep "$i" ); done)
  targIsoseq_cluster=$(for i in ${SUBSET_SAMPLES_NAMES[@]}; do echo $(find $WKD_ROOT/3_refine -name "*flnc.report.csv" | grep "$i" ); done)
  
  # copy and differentiate flnc files
  cd $WKD_ROOT/7b_matched_only/8a_merged_cluster/
  mkdir -p flnctargeted flncwhole
  cp ${globIsoseq_cluster[@]} $WKD_ROOT/7b_matched_only/8a_merged_cluster/flncwhole
  cd $WKD_ROOT/7b_matched_only/8a_merged_cluster/flncwhole; for f in * ; do mv -- "$f" "WholeIso$f" ; done
  mv * ..
  cp ${targIsoseq_cluster[@]} $WKD_ROOT/7b_matched_only/8a_merged_cluster/flnctargeted
  cd $WKD_ROOT/7b_matched_only/8a_merged_cluster/flnctargeted; for f in * ; do mv -- "$f" "TargetedIso$f" ; done
  mv * ../
  rmdir $WKD_ROOT/7b_matched_only/8a_merged_cluster/flncwhole $WKD_ROOT/7b_matched_only/8a_merged_cluster/flnctargeted
  
  source activate nanopore
  demux_merged_cluster.py $WKD_ROOT/7b_matched_only/8a_merged_cluster $WKD_ROOT/7b_matched_only/8a_merged_cluster/$samplename.clustered.cluster_report.csv -o $samplename
  demux_cupcake_collapse.py $WKD_ROOT/7b_matched_only/8c_collapse/$samplename.read_stat.txt $WKD_ROOT/7b_matched_only/8a_merged_cluster/$samplename.clustered.demuxed.cluster_report.csv\
    --dataset=isoseq -o $samplename
}

run_cupcake

# sqanti3
echo "Running SQANTI3..."
mkdir -p $WKD_ROOT/7b_matched_only/8d_sqanti3; cd $WKD_ROOT/7b_matched_only/8d_sqanti3
source activate sqanti2_py3
python $SQANTI3_DIR/sqanti3_qc.py $WKD_ROOT/7b_matched_only/8c_collapse/${samplename}.gff \
  $GENOME_GTF $GENOME_FASTA -t 30 --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA \
  --genename --isoAnnotLite --gff3 $GFF3 --skipORF --report skip --fl_count $WKD_ROOT/7b_matched_only/8c_collapse/${samplename}_fl_count.csv &> ${samplename}_sqanti_qc.log


# sqanti3 filter 
filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_classification.txt \
  --faa=${samplename}_corrected.fasta \
  --gtf=${samplename}_corrected.gtf \
  -j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log
