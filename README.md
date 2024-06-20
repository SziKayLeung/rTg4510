# Transcriptome sequencing of rTg4510 AD mouse model

This directory is a repository of scripts pertaining to whole transcriptome profiling and targeted sequencing of rTg4510 mice, and form the source of the paper: [Long-read transcript sequencing identifies differential isoform expression in the entorhinal cortex in a transgenic model of tau pathology](https://www.biorxiv.org/content/10.1101/2023.09.20.558220v1) by SK.Leung, A.R.Jeffries,..E.Hannon,J.Mill.

## **Summary**
* Following global Iso-Seq transcriptome profiling performed on rTg4510 wild-type (WT) and transgenic (TG) mice using Sequel (PacBio), data from each sample (1 per flow-cell) was merged and processed using Iso-Seq3 pipeline, aligned to reference genome (mm10) and annotated using SQANTI3. 
* Following targeted sequencing of rTg4510 WT and TG mice across 20 dementia- and AD-asociated genes using Sequel (PacBio) and MinION (ONT), respective data was processed using technology-specific pipelines, and merged. 
* Merged transcriptome of 20 AD-genes further characterised using [FICLE](https://github.com/SziKayLeung/FICLE) and long-read proteogenomics pipeline
* Differential expression analysis was performed at a global and targeted level to identify transcript expression changes associated with development of tau pathology in these mice
* Fluorescence-activated nuclei sorting (FANS) followed by ONT nanopore sequencing was used to profile full-length transcripts in purified NeuN+ (neuron-enriched) and NeuN- (glia-enriched) cortical nuclei populations isolated from a subset of mice (n = 8, 2 WT and 2 TG at ages 2 and 8 months).

## **Datasets** 

1. [Whole Iso-Seq Transcriptome (n = 12 samples)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA663877)
2. [Targeted Iso-Seq Transcriptome (n = 24 samples)](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA981131)
3. [Targeted ONT Transcriptome (n = 18 samples)](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA981131)
4. [Whole RNA-Seq Transcriptome (n = 59 samples)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125957)
5. Targeted ONT Transcriptome from FANS dataset (n = 8 samples)

The processed intermediate datafiles and UCSC genome browser tracks (merged mouse targeted data, mouse whole transcriptome Iso-Seq data, mouse sorted data and human AD targeted data) are available for download at [Zenodo](https://zenodo.org/doi/10.5281/zenodo.8101907). 

## Publications
* [Leung et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.09.20.558220v1): Targeted transcriptome sequencing (ONT and PacBio) of rTg4510 mice across 20 dementia-associated genes
* [Leung, Jeffries et al. (2021)](https://www.cell.com/cell-reports/pdf/S2211-1247(21)01504-7.pdf): PacBio Whole Iso-Seq transcriptome sequencing of rTg4510 mice
* [Castanho et al. (2020)](https://www.sciencedirect.com/science/article/pii/S2211124720300887): RNA-Sequencing of rTg4510 mice
