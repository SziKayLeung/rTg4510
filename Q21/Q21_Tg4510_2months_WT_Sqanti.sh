# Q21 - Tg4510 mouse, aged 2months, Entorhinal Cortex, WT
# Date: 10 September 2018 
# Sqanti output 
# After Isoseq3: Collapse_Isoform cupcake script by E.Tseng 
# After GMAP to mouse genome (GRCm38.p4)

# Convert Collapsed/Filtered Isoforms (minimum 2FL passes) from fastq format into fasta format as required for squanti_qc.py  
# https://github.com/lh3/seqtk 
cd /mnt/data1/Szi/Q21
mkdir Sqanti
cd Sqanti
seqtk seq -a /mnt/data1/Szi/Q21/Collapse_Isoform/Q21.collapsed.min_fl_2.filtered.rep.fq > Q21.collapsed.min_fl_2.filtered.fasta

# Q21: Input = Isoform in gtf format
# Reference mouse genome in gtf file 
# Reference mouse genome in fasta file 
# -o = output 
# d = directory for output 

cd /mnt/data1/Szi/Q21/Sqanti

python /mnt/data1/software/PacBio/ConesaLab-sqanti-6927e53e56d2/sqanti_qc.py -g \
/mnt/data1/Szi/Q21/Collapse_Isoform/Q21.collapsed.min_fl_2.filtered.gff \
/mnt/data1/Szi/reference/GRCm38.p4.gtf \
/mnt/data1/Szi/reference/GRCm38.p4.genome.fa \
-o Q21 -d /mnt/data1/Szi/Q21/Sqanti  &> sqanti_qc_outfile.txt

# classification text from sqanti_QC command 
# isoform input in fasta format = original isoform from collapsed (prior to Sqanti_QC)
python /mnt/data1/software/PacBio/ConesaLab-sqanti-6927e53e56d2/sqanti_filter.py \
/mnt/data1/Szi/Q21/Sqanti/Q21_classification.txt \
--isoforms Q21.collapsed.min_fl_2.filtered.fasta \
-d /mnt/data1/Szi/Q21/Sqanti &> sqanti_filter_outfile.txt