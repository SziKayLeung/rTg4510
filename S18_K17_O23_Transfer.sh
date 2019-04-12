#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=2:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

# Date: 12th April 2019 
# Transfer S18, K17, O23 Files from ISCA to ISCA personal scratch drive (to work from as path only recently established)

# L22_Tg4510_TG_8months
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/3_C12/m54082_190306_083150.subreads.bam /gpfs/ts0/scratch/sl693/Isoseq3/rawdata 
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/3_C12/m54082_190306_083150.subreads.bam.pbi /gpfs/ts0/scratch/sl693/Isoseq3/rawdata 
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/3_C12/m54082_190306_083150.subreadset.xml /gpfs/ts0/scratch/sl693/Isoseq3/rawdata 

# K18_Tg4510_WT_2months
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/4_D12/m54082_190307_045507.subreads.bam /gpfs/ts0/scratch/sl693/Isoseq3/rawdata
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/4_D12/m54082_190307_045507.subreads.bam.pbi /gpfs/ts0/scratch/sl693/Isoseq3/rawdata
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190304_153723/4_D12/m54082_190307_045507.subreadset.xml /gpfs/ts0/scratch/sl693/Isoseq3/rawdata

# S18_Tg4510_TG_2mos
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreads.bam /gpfs/ts0/scratch/sl693/Isoseq3/rawdata
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreads.bam.pbi /gpfs/ts0/scratch/sl693/Isoseq3/rawdata
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/2_D11/m54082_190404_101400.subreadset.xml /gpfs/ts0/scratch/sl693/Isoseq3/rawdata 


# K17_Tg4510_WT_2mos
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreads.bam /gpfs/ts0/scratch/sl693/Isoseq3/rawdata
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreads.bam.pbi /gpfs/ts0/scratch/sl693/Isoseq3/rawdata
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190403_131312/3_E11/m54082_190405_063832.subreadset.xml /gpfs/ts0/scratch/sl693/Isoseq3/rawdata


# O23_Tg4510
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreads.bam /gpfs/ts0/scratch/sl693/Isoseq3/rawdata
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreads.bam.pbi /gpfs/ts0/scratch/sl693/Isoseq3/rawdata
cp /bioseqfs/PacBio/Sequel/2801/r54082_20190401_162910/1_A01/m54082_190401_165425.subreadset.xml /gpfs/ts0/scratch/sl693/Isoseq3/rawdata

echo transfer done