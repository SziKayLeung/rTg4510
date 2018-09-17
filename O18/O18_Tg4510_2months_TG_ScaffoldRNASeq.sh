cd /mnt/data1/Szi/O18
mkdir Scaffold
cd Collapse_Isoform_2
# copy both fastq file and gff for index
cp test.collapsed.min_fl_2.filtered.gff ../Scaffold
cp test.collapsed.min_fl_2.filtered.rep.fq ../Scaffold

# convert fastq to fasta file 
cd /mnt/data1/Szi/O18/Scaffold
mv test.collapsed.min_fl_2.filtered.rep.fq ./O18.collapsed.min_fl_2.filtered.rep.fq
seqtk seq -a O18.collapsed.min_fl_2.filtered.rep.fq > O18.collapsed.min_fl_2.filtered.fa

# index pacbio long reads using bowtie2
/mnt/data1/programs/bowtie2-2.3.1/bowtie2-build  O18.collapsed.min_fl_2.filtered.fa O18.collapsed.min_fl_2.filtered &> index.outfile.txt
