##############################################################################################################
# Date: 16th April 2019
# Tabulate ccs output reports into merged final output
#############################################################################################################
import os 
import pandas as pd
from glob import glob
from pathlib import Path

# setting working directory
path = Path(r'/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Isoseq3/LIMA')
os.chdir(path)
# check in working directory
Path.cwd()

# list all files with ccs.report 
filenames = glob('*lima.summary')
filenames

# loop txt.files into one list 
dataframes = [pd.read_csv(f,sep='\t',header = None) for f in filenames]
dataframes

# Generate sample names based on the names of files, extracting the first charachter of the strings delimited by ".""
samples = []
for f in filenames : 
    f = f.split('.')[0]
    samples.append(str(f))
    print(samples)

# To split column into several columns as separated with commas, and save to df1 column
# To label the columns as the sample names (list saved from above)
# Set the index of the dataframe as the descriptions
mod = []
count = 0 
for df in dataframes :
    df.columns = ['Description']
    df1 = df['Description'].str.split(':', expand=True)
    print(len(df1))
    df1.columns = ['Description', samples[count]]
    df1 = df1.set_index('Description')
    mod.append(df1)
    count += 1

# Concentenate dataframes
final = pd.concat(mod, axis = 1)
final.to_csv('LIMA.SUMMARY.csv')
