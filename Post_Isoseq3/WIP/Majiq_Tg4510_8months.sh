#!/bin/bash
# Script to run Majiq build on Tg4510 8months Ctrl vs AD 
# Date: 3rd - 4th December 2018 

## MAJIQ installation
## Created conda environment to install python 3.6 (Majiq Requires >= python 3.5) 
condal create -n py36 
source activate py36 
## MAJIQ installation: Majiq version 1.1.7a 
pip install cython pysam numpy 
conda create -n py36 -c git+https://bitbucket.org/biociphers/majiq_stable.git#egg=majiq

# Build MAJIQ 
# Package installed in: /home/sLeung/.conda/envs/py36/lib/python3.6/site-packages/majiq
cd /home/sLeung/.conda/envs/py36/lib/python3.6/site-packages/majiq
chmod a+x run_majiq.py 
python run_majiq.py

## Majiq Build (version Majiq 1.1.7a)
# Download mouse gff3 reference database 
cd /mnt/data1/Szi/reference 
mkdir MAJIQ
cd MAJIQ 
wget https://majiq.biociphers.org/download/DB/ensembl.mm10.gff3
# Create configuration file # 125bp read length from fastqc 
cd /mnt/data1/Szi
mkdir MAJIQ_Tg4510
cd MAJIQ_Tg4510
nano Tg4510_8mos.ini
##[info]
## readlen=125
## samdir=/mnt/data1/Szi/Mapped_RNA-seq
## genome=mm10
## strandness=None
## [experiments]
## Tg4510_AD_8months=K24Aligned.out.sorted,L22Aligned.out.sorted,M20Aligned.out.sorted,O24Aligned.out.sorted,Q20Aligned.out.sorted,S24Aligned.out.sorted,T22Aligned.out.sorted
## Tg4510_Ctrl_8months=K23Aligned.out.sorted,L21Aligned.out.sorted,M19Aligned.out.sorted,O23Aligned.out.sorted,P21Aligned.out.sorted,Q19Aligned.out.sorted,S23Aligned.out.sorted,T21Aligned.out.sorted

# majiq build <transcript list.gff3> -c <configuration file> -j <number of threads> -o <build outdir>
majiq build /mnt/data1/Szi/reference/MAJIQ/ensembl.mm10.gff3 \
-c /mnt/data1/Szi/MAJIQ_Tg4510/Tg4510_8mos.ini \
-j 32 -o /mnt/data1/Szi/MAJIQ_Tg4510

# majiq weights <files.majiq><replicate1.majiq> -j <number of threads> -n <name of groups> 
cd /mnt/data1/Szi/MAJIQ_Tg4510/
mkdir weights  
cd weights
majiq weights K24Aligned.out.sorted.majiq L22Aligned.out.sorted.majiq M20Aligned.out.sorted.majiq O24Aligned.out.sorted.majiq Q20Aligned.out.sorted.majiq S24Aligned.out.sorted.majiq T22Aligned.out.sorted.majiq \
-j 16 -o /mnt/data1/Szi/MAJIQ_Tg4510/weights -n Tg4510_AD_8months &> weights.output.txt

majiq weights K23Aligned.out.sorted.majiq L21Aligned.out.sorted.majiq M19Aligned.out.sorted.majiq O23Aligned.out.sorted.majiq P21Aligned.out.sorted.majiq Q19Aligned.out.sorted.majiq S23Aligned.out.sorted.majiq T21Aligned.out.sorted.majiq \
-j 16 -o /mnt/data1/Szi/MAJIQ_Tg4510/weights_ctrl -n Tg4510_Ctrl_8months &> weights_ctrl.output.txt

# majiq deltapsi -grp1 <files.majiq> -grp2 <files.majiq> -j <number of threads> -o <build outdir> 
cd /mnt/data1/Szi/MAJIQ_Tg4510
majiq deltapsi -grp1 K24Aligned.out.sorted.majiq L22Aligned.out.sorted.majiq M20Aligned.out.sorted.majiq O24Aligned.out.sorted.majiq Q20Aligned.out.sorted.majiq S24Aligned.out.sorted.majiq T22Aligned.out.sorted.majiq \
-grp2 K23Aligned.out.sorted.majiq L21Aligned.out.sorted.majiq M19Aligned.out.sorted.majiq O23Aligned.out.sorted.majiq P21Aligned.out.sorted.majiq Q19Aligned.out.sorted.majiq \
-j 16 -o /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi -n Tg4510_AD_8months Tg4510_Ctrl_8months &> delta.psi.outputs.txt

# voila deltapsi <dpsi outdir>/<cond1_id>_<cond2_id>.deltapsi.voila -s <build outdir>/splicegraph.sql -o <voila outdir>
voila deltapsi /mnt/data1/Szi/MAJIQ_Tg4510/delta_psi/Tg4510_AD_8months_Tg4510_Ctrl_8months.deltapsi.voila -s /mnt/data1/Szi/MAJIQ_Tg4510/splicegraph.sql -o /mnt/data1/Szi/MAJIQ_Tg4510/voila &> voila.psi.outputs.txt

