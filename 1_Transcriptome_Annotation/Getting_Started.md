## Getting Started

1.  Create/load conda environment to run the pipeline in ISCA, which contain functions that call python scripts with dependencies. All installation instructions found in [CreatingCondaEnv_Installation.sh](https://git.exeter.ac.uk/sl693/general/blob/master/Reference/CreatingCondaEnv_Installation.sh).

```
module load Miniconda2
source activate isoseq3     # all necessary dependencies to run IsoSeq3
source activate sqanti2_py3 # all necessary dependencies to run RNASeq and Post-IsoSeq (Cupcake and SQANTI2) in python3
source activate sqanti3     # Post-Isoseq (Cupcake and SQANTI2) in python3

```

2. Download all the genome.gtf file from GENCODE and build an index using STAR, GMAP and Minimap2 for alignment downstream
3. Create a fasta sequence for the cDNA primer sequences for demultiplexing downstream
4. Download mm9 FANTOM5 CAGE peak and use a liftover to mm10 for input into SQANTI2 ([Installation Instructions](https://git.exeter.ac.uk/sl693/general/blob/master/Reference/CreatingCondaEnv_Installation.sh)) 
