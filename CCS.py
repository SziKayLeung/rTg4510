import os 
import pandas as pd
from glob import glob
from pathlib import Path

# setting working directory
path = Path(r'C:\Users\Kay\Documents\CCS')
os.chdir(path)
# check in working directory
pathlib.Path.cwd()

filenames = glob('*ccs_report.txt')
filenames

dataframe = [pd.read_csv(f,sep='\t') for f in filenames]
dataframe

len(dataframe.columns)

samples = []
for f in filenames : 
    f = f.split('_')[0]
    samples.append(str(f))
    print(samples)

mod = []
count = 0 
for df in dataframe :
    df1 = df['ZMW Yield'].str.split(',', expand=True)
    df1.columns = ['ZMW Yield', samples[count], samples [count]]
    df1 = df1.set_index('ZMW Yield')
    mod.append(df1)
    count += 1

df_col = pd.concat([mod[0],mod[1],mod[2]],axis = 1)