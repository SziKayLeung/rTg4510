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
#SBATCH --error=4_blastvsBDR.e
#SBATCH --output=4_blastvsBDR.o


# 05/02/2024: Blast against BDR targeted dataset (Sqanti3)

##-------------------------------------------------------------------------

#Blast_seq2seq <working directory> <sample.fasta.for.db> <output_db_name> <sample.fasta.for.blast> <blast.output.name> <build> <threshold>
Blast_seq2seq(){
  
  ### Prepared Reference Genome for BLAST
  module load Miniconda2
  source activate nanopore
  cd $1
  if [[ $6 == "build" ]]; then 
  echo "Create blast database with fasta sequence: $2"
  makeblastdb -in $2 -dbtype nucl -out $3
  echo "Output file successfully generated: $2"
  
  ### Blast Probes to Reference Genome
  # https://angus.readthedocs.io/en/2016/running-command-line-blast.html
  # http://envgen.nox.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf
  echo "Blast Probes to indexed database with fasta sequence: $3"
  if [[ $7 == "high" ]]; then
  blastn -query $4 -db $3 -out $5 -outfmt 6 -evalue 1e-5 
  #  -outfmt "7 qacc sacc evalue qstart qend sstart send" -evalue 1e-5  
  elif [[ $7 == "low" ]]; then
  echo "Using low threshold for evalue"
  blastn -query $4 -db $3 -out $5 -outfmt 6 -evalue 1e-2
  else
    echo "Threshold required as high or low"
  fi
  
  
  # column headers:https://molevol.mbl.edu/index.php/BLAST_UNIX_Tutorial
  #echo -e "Query_ID\Subject_ID\%_Identity\alignment_length\mismatches\gap\q.start\q.end\s.start\s.end\e_value\bit_score" | cat - $5 > Final_$5
  else 
    echo "Blast Probes to indexed database with fasta sequence: $3"
  if [[ $7 == "high" ]]; then
  blastn -query $4 -db $3 -out $5 -outfmt 6 -evalue 1e-5 
  #  -outfmt "7 qacc sacc evalue qstart qend sstart send" -evalue 1e-5  
  elif [[ $7 == "low" ]]; then
  echo "Using low threshold for evalue"
  blastn -query $4 -db $3 -out $5 -outfmt 6 -evalue 10
  else
    echo "Threshold required as high or low"
  fi
  fi  
  
}

outputDir=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/4_characterise/BDRBlast
mousefasta=/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/G_Merged_Targeted/2_sqanti3/all_iso_ont_collapsed.filtered_counts_filtered.fa
humanfasta=/lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/5_cupcake/7_sqanti3/ontBDR_collapsed.filtered_counts_filtered.fa
Blast_seq2seq ${outputDir} ${mousefasta} allIsoOnt.db ${humanfasta} allIsoOntBDRHigh.txt build high 
Blast_seq2seq ${outputDir} ${mousefasta} allIsoOnt.db ${humanfasta} allIsoOntBDRLow.txt notbuild low 

