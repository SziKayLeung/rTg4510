#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Transfer_Rawdata_Batch2_3_Isoseq.out

# 22/04/2020: Transfer Targeted Sequencing (Batch 1) Raw data from Project_10040 sequencing service to ISCA 
# 26/10/2020: Transfer Targeted Sequencing (Batch 2, Batch 3) Raw data (Iso-Seq) from Project 10175 sequencing service to ISCA

# To transfer new file 
# 1. Include in "project_name" list the name of new file 
# 2. Create new variable for ftp link of new file 
# 3. Increase counter in the until loop i.e. -gt <number_of_projects>

####################### Define variables 
# Project Name
# project_name=(<project_name_1> <project_name_2>)
project_name=(Project_10040 Project_10175)

# Directory Paths 
Targeted_rawdata=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Targeted/rawdata

# ftp link of raw sequencing data 
# !!!Ensure ftp link order refer to project_name order!!!
# FTP=(<Project_name_raw_1> <Project_name_raw_2>)
Project_10040_raw="ftp://Project_10040:B9jeZ3TE49j2@ftp1.sequencing.exeter.ac.uk/"
Project_10175_raw="ftp://Project_10175:nFreXBFceJtl@ftp1.sequencing.exeter.ac.uk/"
FTP=($Project_10040_raw $Project_10175_raw)


####################### Apply Loop

project_num=0 # counter
until [ $project_num -gt 2 ]; do
    echo project_num: $project_num
    
    working_project_name=${project_name[project_num]}
    working_project_dir=$Targeted_rawdata/${project_name[project_num]}
    working_ftp=${FTP[project_num]}
    
    echo "Processing: $working_project_name"

    if [ -d $working_project_dir ]; then 
        ### Take action if $working_project_dir exists ###
        echo "Already transferred $working_project_name to $working_project_dir"
    else
         ###  Control will jump here if $working_project_dir does NOT exists ###
        echo "Transfer $working_ftp to $working_project_dir"
        
        cd $Targeted_rawdata
        
        if wget -m $working_ftp &> $working_project_name"_transfer_log" ; then  
            echo "Extract $working_project_name sucessful"
        fi
        
        mv ftp1.sequencing.exeter.ac.uk/ $working_project_name
        echo "Log file of $working_project_name"_transfer_log""
        head $working_project_name"_transfer_log"
        tail $working_project_name"_transfer_log"
    fi
  
  ((project_num++))
done
