# Q21 - Tg4510 mouse, aged 2months, Entorhinal Cortex, WT
# Date: 06 September 2018 
# After Isoseq3: Collapse_Isoform cupcake script by E.Tseng 
# After GMAP to mouse genome (GRCm38.p4)


# Collapse Isoform after mapping to GMAP 
# remove reads not mapped to GMAP, have low coverage, low identity  
cd /mnt/data1/Szi/Q21  
mkdir Collapse_Isoform 
cd Collapse_Isoform

# 1. To create cluster report for abundance.txt (cluster report = list of the clustered full length and non-full length)
# /home/sLeung/cDNA_Cupcake/post_isoseq_cluster/isoseq3_make_cluster_report.py = path to script (Note updated from E.Tseng September 2018 after formatting ID)  
# /mnt/data1/Szi/Q21/Isoseq3/Q21.polished.bam = path to polished.bam (Isoseq3 output) 
python /home/sLeung/cDNA_Cupcake/post_isoseq_cluster/isoseq3_make_cluster_report.py \
/mnt/data1/Szi/Q21/Isoseq3/Q21.polished.bam &> Q21_cluster_logfile.txt


# 2. Collapse Isoforms
# /home/sLeung/cDNA_Cupcake/cupcake/tofu/collapse_isoforms_by_sam.py = path to script 
# --input /mnt/data1/Szi/Q21/Isoseq3/Q21.polished.hq.fastq --fq = input fastq after polishing 
# -s /mnt/data1/Szi/Q21/GMAP/hq_isoforms.fastq.sorted.sam = sorted sam file from GMAP 
# -o Q21 = output (for ease sake)

python /home/sLeung/cDNA_Cupcake/cupcake/tofu/collapse_isoforms_by_sam.py \
--input /mnt/data1/Szi/Q21/Isoseq3/Q21.polished.hq.fastq --fq \
-s /mnt/data1/Szi/Q21/GMAP/hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o Q21 &>  Q21_collapse_logfile.txt


# 3. Create Abundance Script of full-length transcripts  
cd /mnt/data1/Szi/Q21/Collapse_Isoform
# /home/sLeung/cDNA_Cupcake/cupcake/tofu/get_abundance_post_collapse.py = path to script 
# Q21.collapsed = name of collapsed files created from Collapse_Isoform (note output from before i.e. Q21)
# cluster_report.csv = cluster report created from line 13 
python /home/sLeung/cDNA_Cupcake/cupcake/tofu/get_abundance_post_collapse.py Q21.collapsed cluster_report.csv &>  Q21_abundance_logfile.txt


# 4. Filter by counts of 2 full length passes 
cd /mnt/data1/Szi/Q21/Collapse_Isoform 
python /home/sLeung/cDNA_Cupcake/cupcake/tofu/filter_by_count.py /mnt/data1/Szi/Q21/Collapse_Isoform/Q21.collapsed --min_count=2 --dun_use_group_count &> Q21.filter.txt


# 5. remove degraded isoforms (default setting) 
cd /mnt/data1/Szi/Q21/Collapse_Isoform 
python /home/sLeung/cDNA_Cupcake/cupcake/tofu/filter_away_subset.py Q21.collapsed.min_fl_2 &>  Q21_filtered_logfile.txt