#!/bin/bash

## MAJIQ installation
## Created conda environment to install python 3.6 (Majiq Requires >= python 3.5) 
#condal create -n py36 
#source activate py36 
## MAJIQ installation: Majiq version 1.1.7a 
#pip install cython pysam numpy 
#conda create -n py36 -c git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq

# Build MAJIQ 
# Package installed in: /home/sLeung/.conda/envs/py36/lib/python3.6/site-packages/majiq
#cd /home/sLeung/.conda/envs/py36/lib/python3.6/site-packages/majiq
#chmod a+x run_majiq.py 
#python run_majiq.py

## Majiq Build (version Majiq 1.1.7a)
# Download mouse gff3 reference database 
#cd /mnt/data1/Szi/reference 
#mkdir MAJIQ
#cd MAJIQ 
#wget https://majiq.biociphers.org/download/DB/ensembl.mm10.gff3

## Create configuration file # 125bp read length from fastqc 
## cd /mnt/data1/Szi/MAJIQ_Tg4510
## nano Tg4510.ini
##[info]
##readlen=125
##samdir=/mnt/data1/Szi/Mapped_RNA-seq
##genome=mm10
##strandness=None
##[experiments]
##Tg4510_Ctrl_2months=K17Aligned.out.sorted,L23Aligned.out.sorted,M21Aligned.out.sorted,P23Aligned.out.sorted,Q21Aligned.out.sorted,S17Aligned.out.sorted,T23Aligned.out.sorted
##Tg4510_AD_2months=K18Aligned.out.sorted,L24Aligned.out.sorted,M22Aligned.out.sorted,O18Aligned.out.sorted,P24Aligned.out.sorted,Q22Aligned.out.sorted,S18Aligned.out.sorted,T24Aligned.out.sorted
##Tg4510_Ctrl_4months=K19Aligned.out.sorted,L17Aligned.out.sorted,M23Aligned.out.sorted,O19Aligned.out.sorted,P17Aligned.out.sorted,Q23Aligned.out.sorted,S19Aligned.out.sorted,T17Aligned.out.sorted
##Tg4510_AD_4months=K20Aligned.out.sorted,L18Aligned.out.sorted,M24Aligned.out.sorted,O20Aligned.out.sorted,P18Aligned.out.sorted,Q24Aligned.out.sorted,S20Aligned.out.sorted,T18Aligned.out.sorted
##Tg4510_Ctrl_6months=K21Aligned.out.sorted,M17Aligned.out.sorted,O21Aligned.out.sorted,P19Aligned.out.sorted,Q17Aligned.out.sorted,S21Aligned.out.sorted,T19Aligned.out.sorted
##Tg4510_AD_6months=M18Aligned.out.sorted,O22Aligned.out.sorted,P20Aligned.out.sorted,Q18Aligned.out.sorted,S22Aligned.out.sorted,T20Aligned.out.sorted
##Tg4510_AD_8months=K24Aligned.out.sorted,L22Aligned.out.sorted,M20Aligned.out.sorted,O24Aligned.out.sorted,Q20Aligned.out.sorted,S24Aligned.out.sorted,T22Aligned.out.sorted
##Tg4510_Ctrl_8months=K23Aligned.out.sorted,L21Aligned.out.sorted,M19Aligned.out.sorted,O23Aligned.out.sorted,P21Aligned.out.sorted,Q19Aligned.out.sorted,S23Aligned.out.sorted,T21Aligned.out.sorted

## Majiq build 
# majiq build /mnt/data1/Szi/reference/MAJIQ/ensembl.mm10.gff3 -c /mnt/data1/Szi/MAJIQ_Tg4510/Tg4510.ini -j 16 -o /mnt/data1/Szi/MAJIQ_Tg4510 &> Tg4510.build.outputs.txt

cd /mnt/data1/Szi/MAJIQ_Tg4510

# AD_Ctrl_8mos
## majiq deltapsi -grp1 K24Aligned.out.sorted.majiq L22Aligned.out.sorted.majiq M20Aligned.out.sorted.majiq O24Aligned.out.sorted.majiq Q20Aligned.out.sorted.majiq S24Aligned.out.sorted.majiq T22Aligned.out.sorted.majiq \
## -grp2 K23Aligned.out.sorted.majiq L21Aligned.out.sorted.majiq M19Aligned.out.sorted.majiq O23Aligned.out.sorted.majiq P21Aligned.out.sorted.majiq Q19Aligned.out.sorted.majiq \
## -j 16 -o /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/AD_Ctrl_8mos -n Tg4510_AD_8months Tg4510_Ctrl_8months &> delta_psi_AD_Ctrl_8mos.outputs.txt

# AD_Ctrl_2mos
##Tg4510_AD_2months=K18Aligned.out.sorted,L24Aligned.out.sorted,M22Aligned.out.sorted,O18Aligned.out.sorted,P24Aligned.out.sorted,Q22Aligned.out.sorted,S18Aligned.out.sorted,T24Aligned.out.sorted
##Tg4510_Ctrl_2months=K17Aligned.out.sorted,L23Aligned.out.sorted,M21Aligned.out.sorted,P23Aligned.out.sorted,Q21Aligned.out.sorted,S17Aligned.out.sorted,T23Aligned.out.sorted
## majiq deltapsi -grp1 K18Aligned.out.sorted.majiq L24Aligned.out.sorted.majiq M22Aligned.out.sorted.majiq O18Aligned.out.sorted.majiq P24Aligned.out.sorted.majiq Q22Aligned.out.sorted.majiq S18Aligned.out.sorted.majiq T24Aligned.out.sorted.majiq \
## -grp2 K17Aligned.out.sorted.majiq L23Aligned.out.sorted.majiq M21Aligned.out.sorted.majiq P23Aligned.out.sorted.majiq Q21Aligned.out.sorted.majiq S17Aligned.out.sorted.majiq T23Aligned.out.sorted.majiq \
## -j 16 -o /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/AD_Ctrl_2mos -n Tg4510_AD_2months Tg4510_Ctrl_2months &> delta_psi_AD_Ctrl_2mos.outputs.txt

# AD_Ctrl_4mos
##Tg4510_AD_4months=K20Aligned.out.sorted,L18Aligned.out.sorted,M24Aligned.out.sorted,O20Aligned.out.sorted,P18Aligned.out.sorted,Q24Aligned.out.sorted,S20Aligned.out.sorted,T18Aligned.out.sorted
##Tg4510_Ctrl_4months=K19Aligned.out.sorted,L17Aligned.out.sorted,M23Aligned.out.sorted,O19Aligned.out.sorted,P17Aligned.out.sorted,Q23Aligned.out.sorted,S19Aligned.out.sorted,T17Aligned.out.sorted
## majiq deltapsi -grp1 K20Aligned.out.sorted.majiq L18Aligned.out.sorted.majiq M24Aligned.out.sorted.majiq O20Aligned.out.sorted.majiq P18Aligned.out.sorted.majiq Q24Aligned.out.sorted.majiq S20Aligned.out.sorted.majiq T18Aligned.out.sorted.majiq \
## -grp2 K19Aligned.out.sorted.majiq L17Aligned.out.sorted.majiq M23Aligned.out.sorted.majiq O19Aligned.out.sorted.majiq P17Aligned.out.sorted.majiq Q23Aligned.out.sorted.majiq S19Aligned.out.sorted.majiq T17Aligned.out.sorted.majiq \
## -j 16 -o /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/AD_Ctrl_4mos -n Tg4510_AD_4months Tg4510_Ctrl_4months &> delta_psi_AD_Ctrl_4mos.outputs.txt

# AD_Ctrl_6mos
##Tg4510_AD_6months=M18Aligned.out.sorted,O22Aligned.out.sorted,P20Aligned.out.sorted,Q18Aligned.out.sorted,S22Aligned.out.sorted,T20Aligned.out.sorted
##Tg4510_Ctrl_6months=K21Aligned.out.sorted,M17Aligned.out.sorted,O21Aligned.out.sorted,P19Aligned.out.sorted,Q17Aligned.out.sorted,S21Aligned.out.sorted,T19Aligned.out.sorted
## majiq deltapsi -grp1 M18Aligned.out.sorted.majiq O22Aligned.out.sorted.majiq P20Aligned.out.sorted.majiq Q18Aligned.out.sorted.majiq S22Aligned.out.sorted.majiq T20Aligned.out.sorted.majiq \
## -grp2 K21Aligned.out.sorted.majiq M17Aligned.out.sorted.majiq O21Aligned.out.sorted.majiq P19Aligned.out.sorted.majiq Q17Aligned.out.sorted.majiq S21Aligned.out.sorted.majiq T19Aligned.out.sorted.majiq \
## -j 16 -o /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/AD_Ctrl_6mos -n Tg4510_AD_6months Tg4510_Ctrl_6months &> delta_psi_AD_Ctrl_6mos.outputs.txt

voila deltapsi /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/AD_Ctrl_8mos/Tg4510_AD_8months_Tg4510_Ctrl_8months.deltapsi.voila -s /mnt/data1/Szi/MAJIQ_Tg4510/splicegraph.sql \
-o /mnt/data1/Szi/MAJIQ_Tg4510/voila/voila_8mos &> voila.AD_Ctrl_8mos.outputs.txt

voila deltapsi /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/AD_Ctrl_4mos/Tg4510_AD_4months_Tg4510_Ctrl_4months.deltapsi.voila -s /mnt/data1/Szi/MAJIQ_Tg4510/splicegraph.sql \
-o /mnt/data1/Szi/MAJIQ_Tg4510/voila/voila_2mos &> voila.AD_Ctrl_4mos.outputs.txt

voila deltapsi /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/AD_Ctrl_6mos/Tg4510_AD_6months_Tg4510_Ctrl_6months.deltapsi.voila -s /mnt/data1/Szi/MAJIQ_Tg4510/splicegraph.sql \
-o /mnt/data1/Szi/MAJIQ_Tg4510/voila/voila_6mos &> voila.AD_Ctrl_6mos.outputs.txt

voila deltapsi /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/AD_Ctrl_2mos/Tg4510_AD_2months_Tg4510_Ctrl_2months.deltapsi.voila -s /mnt/data1/Szi/MAJIQ_Tg4510/splicegraph.sql \
-o /mnt/data1/Szi/MAJIQ_Tg4510/voila/voila_2mos &> voila.AD_Ctrl_2mos.outputs.txt