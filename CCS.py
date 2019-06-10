#!/usr/bin/env python3
#PBS -V # export all environment variables to the batch job.
#PBS -q mrchq # submit to the serial queue
#PBS -l walltime=1:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=32 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion

##############################################################################################################
# Date: 16th April 2019
# Tabulate ccs output reports into merged final output
#############################################################################################################

import os 
import pandas as pd
from glob import glob
from pathlib import Path

# setting working directory
path = Path(r'/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/WT8_ISOSEQ/IsoSeq3.1.2/CCS')
os.chdir(path)
# check in working directory
Path.cwd()

# list all files with ccs.report 
filenames = glob('*ccs_report.txt')
filenames

# loop txt.files into one list 
dataframes = [pd.read_csv(f,sep='\t') for f in filenames]
dataframes

# Generate sample names based on the names of files, extracting the first charachter of the strings delimited by ".""
samples = []
for f in filenames : 
    f = f.split('_')[0]
    print(f)
    samples.append(str(f))
    print(samples)

# To split column into several columns as separated with commas, and save to df1 column
# To label the columns as the sample names (list saved from above)
# Set the index of the dataframe as the descriptions
mod = []
count = 0 
for df in dataframes :
    df1 = df['ZMW Yield'].str.split(',', expand=True)
    df1.columns = ['ZMW Yield', samples[count], samples [count]]
    df1 = df1.set_index('ZMW Yield')
    mod.append(df1)
    count += 1

# Concentenate dataframes
final = pd.concat([mod[0],mod[1],mod[2],mod[3],mod[4],mod[5],mod[6],mod[7]],axis = 1)
final.to_csv("CCS_output.csv")