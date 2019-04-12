#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 12th April 2019
# run GMAP for all Tg4510 samples: L22, K18, O23, S18, K17 
# already Isoseq3
#############################################################################################################
module load cDNA_Cupcake/6.9-intel-2017b-Python-2.7.14

#############################################################################################################
# directory for saving output
cd /gpfs/ts0/scratch/sl693/WholeTranscriptome/ToFU
# determine path directory for input data
GMAP=/gpfs/ts0/scratch/sl693/WholeTranscriptome/Isoseq3/GMAP

SAMPLES_NAMES=(L22 K18 S18 K17 O23)
#############################################################################################################
# Collapse
for sample in "${SAMPLES_NAMES[@]}"; do 
  echo "Processing $sample file for collapse"
  collapse_isoforms_by_sam.py --input $GMAP/"$sample.polished.hq.fastq" --fq -s $GMAP/$sample.hq_isoforms.fastq.sorted.sam \
  --dun-merge-5-shorter -o $sample   
  echo "Collapse $sample successful"
done

# Create Abundance Script of full-length transcripts  


python /home/sLeung/cDNA_Cupcake/cupcake/tofu/get_abundance_post_collapse.py MouseCTX.collapsed /mnt/data1/Szi/Mouse_CTX_Isoseq3.1/polish/polished.cluster_report.csv

# Filter by counts of 2 full length passes 
python /home/sLeung/cDNA_Cupcake/cupcake/tofu/filter_by_count.py MouseCTX.collapsed --min_count=2 --dun_use_group_count 

# Remove degraded isoforms (default setting) 
python /home/sLeung/cDNA_Cupcake/cupcake/tofu/filter_away_subset.py MouseCTX.collapsed.min_fl_2

conda deactivate

# Sqanti

cd /mnt/data1/Szi/Mouse_CTX_Isoseq3.1/Mapped_GRCm38.p4
mkdir SQANTI

python /mnt/data1/software/PacBio/ConesaLab-sqanti-6927e53e56d2/sqanti_qc.py -g \
/mnt/data1/Szi/Mouse_CTX_Isoseq3.1/Mapped_GRCm38.p4/MouseCTX.collapsed.min_fl_2.filtered.gff \
/mnt/data1/Szi/reference/GRCm38.p4.gtf \
/mnt/data1/Szi/reference/GRCm38.p4.genome.fa \
-o MouseCTX -d /mnt/data1/Szi/Mouse_CTX_Isoseq3.1/Mapped_GRCm38.p4/SQANTI
