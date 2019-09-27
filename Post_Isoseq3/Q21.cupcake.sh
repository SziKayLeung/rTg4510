cd /mnt/data1/Szi/Q21

gmap -D /mnt/data1/Szi/reference -d GRCm38.p4 -f samse -n 0 -t 12 \
   -z sense_force /mnt/data1/Szi/Q21/Q21.polished.hq.fastq\
> hq_isoforms.fastq.sam

sort -k 3,3 -k 4,4n hq_isoforms.fastq.sam > hq_isoforms.fastq.sorted.sam

mkdir Collapse_Isoform
cd Collapse_Isoform

collapse_isoforms_by_sam.py --input /mnt/data1/Szi/Q21/Q21.polished.hq.fastq --fq \
   -s /mnt/data1/Szi/Q21/hq_isoforms.fastq.sorted.sam  --dun-merge-5-shorter -o Q21
