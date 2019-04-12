cd /mnt/data1/Szi/O18
mkdir O18_RNAseq_align
cd O18_RNAseq_align
STAR --runThreadN 32 --genomeDir /mnt/data1/Szi/reference/STAR \
--readFilesIn /mnt/data1/isabel/RNA-seq/analysis_isabel_new/Tg4510_filtered/O18/O18_S36_R1_001.fastq.filtered \
/mnt/data1/isabel/RNA-seq/analysis_isabel_new/Tg4510_filtered/O18/O18_S36_R2_001.fastq.filtered

echo O18 STAR Alignment Successful

samtools view -S -b Aligned.out.sam > Aligned.out.bam
samtools sort Aligned.out.bam Aligned.out.sorted
echo SAMTOOLS sort Successful
samtools index Aligned.out.sorted.bam
echo SAMTOOLS index Successful
samtools flagstat Aligned.out.sorted.bam > mappingstats.txt

mv Aligned.out.sorted.bam O18_Tg4510_2mos_TG.sorted.bam
mv Aligned.out.sorted.bam.bai O18_Tg4510_2mos_TG.sorted.bam.bai
rm Aligned.out.sam
