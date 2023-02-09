## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##



## ---------- Source function and config files -----------------

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/rTg4510/Paper_Figures/"
source(paste0(SC_ROOT, "0_source_functions.R"))
source(paste0(SC_ROOT, "rTg4510_config.R"))
source(paste0(SC_ROOT,"bin/draw_heatmap_gene_level.R"))


# heatmap 
draw_heatmap_gene_level(loaded$glob_iso$results_gene, annotated$glob_iso$GeneExp,"glob_isoseq",diff="yes")
draw_heatmap_gene_level(loaded$glob_iso$results_gene, annotated$glob_iso$GeneExp,"glob_isoseq",diff="no")


common_diff_genes <- intersect(c(diff_models$glob_iso_siggenes$models$`Model 1 Genotype`$...1,
                                 diff_models$glob_iso_siggenes$models$`Model 2 Genotype+Age`$...1),
                               rnaseq_results$GenotypeDEG$Gene)


dat <- annotated$glob_iso$GeneExp[annotated$glob_iso$GeneExp$associated_gene %in% common_diff_genes,] 

dat <- dat %>% filter(time == 8)
dat <- aggregate(Exp ~ group + associated_gene, data = dat, mean) %>%
  spread(., group, Exp) %>% 
  mutate(log2FCIso = log2(CASE/CONTROL), 
         dataset = "IsoSeq") %>% 
  select(associated_gene, log2FCIso, dataset)

rnaseq_results$GenotypeDEG %>% filter(Gene %in% common_diff_genes) %>% 
  mutate(WT_Mean = word(Mean.WT..SD.,c(1), sep = fixed(" ")),
         TG_Mean = word(Mean.TG..SD.,c(1), sep = fixed(" ")))
