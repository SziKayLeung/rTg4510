#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Aim: selection of samples for IHC experiment (n = 8) by plotting the number of reads detected for the most abundant Trem2 isoform
## --------------------------------


library("ggplot2")
library("dplyr")
library("ggrepel")
library("cowplot")

# full-length read counts of PB.20818.54  (all_iso_ont_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered)
trem2 = read.table("/lustre/projects/Research_Project-MRC148213/lsl693/scripts/rTg4510/trem2_reads.txt",header=T)

ONT = trem2 %>% filter(seq == "ONT") %>% 
  ggplot(., aes(x = as.factor(age), y = reads, group = age, label = sample)) + geom_point() + 
  facet_grid(~ phenotype, space = "free", scales = "free") +
  geom_text_repel() + labs(x = "Age (months)", y = "Number of FL reads", title = "ONT") + theme_classic()

Iso = trem2 %>% filter(seq == "Iso.Seq") %>% 
  ggplot(., aes(x = as.factor(age), y = reads, group = age, label = sample)) + geom_point() + 
  facet_grid(~ phenotype, space = "free", scales = "free") +
  geom_text_repel() + labs(x = "Age (months)", y = "Number of FL reads", title = "Iso-Seq") + theme_classic()


plot_grid(ONT, Iso)
