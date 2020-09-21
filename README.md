# Isoseq3_Tg4510

This directory stores scripts pertaining to whole transcriptome Iso-Seq analysis of Tg4510 mice

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

---
## Bioinformatics Pipeline 
1. **Iso-Seq3 QC**: Contain scripts for testing paramaters and for general QC of Iso-Seq3 output. 
    + Determining parameters for number of passes and minimum base accuracy (RQ) for CCS generation from the [number of output reads](https://git.exeter.ac.uk/sl693/IsoSeq3_Tg4510/-/blob/master/Isoseq3/Isoseq3_QC/IsoSeq_Summary_Testing.Rmd) and [ERCC detection](https://git.exeter.ac.uk/sl693/IsoSeq3_Tg4510/-/blob/master/Isoseq3/Isoseq3_QC/IsoSeq_ERCC_Testing.Rmd)
    + [Yield](https://git.exeter.ac.uk/sl693/IsoSeq3_Tg4510/-/blob/master/Isoseq3/Isoseq3_QC/Tg4510_RunStats.Rmd) from Sequel Run output across Tg4510 samples 
    + [Validated](https://git.exeter.ac.uk/sl693/IsoSeq3_Tg4510/-/blob/master/Isoseq3/Isoseq3_QC/Tg4510_hMAPT.R) WT samples did not have any human MAPT (i.e. check samples are wild type)

<br>

2. **Iso-Seq3**: Run Isoseq3 bioinformatics pipeline on [all Tg4510 samples](https://git.exeter.ac.uk/sl693/IsoSeq3_Tg4510/-/blob/master/Isoseq3/Isoseq3_all.sh), [WT Tg4510 samples only](https://git.exeter.ac.uk/sl693/IsoSeq3_Tg4510/-/blob/master/Isoseq3/WT_Whole_Isoseq3.sh) and [TG Tg4510 samples only](https://git.exeter.ac.uk/sl693/IsoSeq3_Tg4510/-/blob/master/Isoseq3/TG_Whole_Isoseq3.sh), from sourcing [function script](https://git.exeter.ac.uk/sl693/general/-/blob/master/IsoSeq/Isoseq3.2.2_Functions.sh)

<br>

3. **Post Iso-Seq3** 
   
---

