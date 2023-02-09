
# check samples where DMP is upregulated
p1 <- RRBS_completebetas %>% filter(position == "chr18:32431597") %>% 
  reshape2::melt() %>% 
  left_join(., RRBS_Phenotype[,c("Sample_ID", "Genotype", "Age_months")], by = c("variable" = "Sample_ID")) %>% 
  filter(Age_months == "8") %>% 
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  ggplot(., aes(x = Genotype, y = value)) + geom_boxplot() + geom_jitter() + 
  labs(x = "Genotype", y = "Methylation") + 
  geom_label_repel(aes(label = variable), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') 


# check samples where isoform is upregulated
p2 <- annotated$ont$Norm_transcounts %>% filter(isoform == "ENSMUST00000234496.1") %>% 
  filter(time == 8) %>%
  mutate(Genotype = factor(ifelse(group == "CONTROL", "WT", "TG"), levels = c("WT","TG"))) %>% 
  ggplot(., aes(x = Genotype, y = value)) + geom_boxplot() + geom_jitter() +
  labs(x = "Genotype", y = "Expression (normalised counts)") + 
  geom_label_repel(aes(label = sample), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') 

plot_grid(p1,p2)


# 
bin1_char.files.names <- list(
  co = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/4_characterise/TargetGenes/Bin1/Stats/Bin1flattened_gencode.csv",
  ES = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/4_characterise/TargetGenes/Bin1/Stats/Bin1_Exonskipping_tab.csv",
  ES_general = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/4_characterise/TargetGenes/Bin1/Stats/Bin1_Exonskipping_generaltab.csv",
  novelexon_co = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/4_characterise/TargetGenes/Bin1/Stats/Bin1_novelexon_coordinates.csv",
  noveexon = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/rTg4510/Merged_Targeted/4_characterise/TargetGenes/Bin1/Stats/Bin1_NovelExon_generaltab.csv"
) 
bin1_char <- lapply(bin1_char.files.names, function(x) read.csv(x))

bin1_char$ES %>% filter(transcript_id == "PB.3915.33_ENSMUST00000234496.1")
# exon 17 
bin1_char$co %>% filter(updated_exon_number == 17)
iso_E17_skipped <- bin1_char$ES %>% filter(ES == "Gencode_17") %>% 
  mutate(ONT_isoform = word(transcript_id, c(2), sep = fixed("_"))) %>% 
  mutate(ONT_isoform = ifelse(grepl("_", transcript_id), word(transcript_id, c(2), sep = fixed("_")),as.character(transcript_id)),
         ONT_isoform = ifelse(grepl("PB",ONT_isoform),"NA",ONT_isoform))
cat("Number of isoforms with also exon 17 skipped:", nrow(iso_E17_skipped))

# 9 isoforms that are also skipped with exon 17 that is upregulated
diff_exon17_skipped <- loaded$ont$results_trans %>% tibble::rownames_to_column(., var = "isoform") %>% filter(isoform %in% iso_E17_skipped$ONT_isoform) 

annotated$ont$Norm_transcounts %>% 
  filter(isoform %in% diff_exon17_skipped$isoform) %>% 
  mutate(Genotype = factor(ifelse(group == "CONTROL", "WT", "TG"), levels = c("WT","TG"))) %>% 
  filter(time == 8) %>%
  ggplot(., aes(x = Genotype, y = log10(value))) + geom_boxplot() + geom_jitter() + 
  facet_grid(~isoform)  + 
  labs(x = "Genotype", y = "Expression (log10 normalised counts)")


# check that the primers unique to ENSMUST00000234496.1 (span across 12 - 18)
# exon 12: 32424820 32424967
# exon 18: 32431983 32432093 
# ENSMUST00000234496.1 contain Gencode 12, and 18, but skipped 13, 14, 15, 16, 17

primer_specific <- bin1_char$ES_general %>% filter(Gencode_12 == "No" & Gencode_18 == "No" & 
                                  Gencode_13 == "Yes" & Gencode_14 == "Yes" &
                                  Gencode_15 == "Yes" & Gencode_16 == "Yes" & Gencode_17 == "Yes") %>% 
  mutate(ONT_isoform = ifelse(grepl("_", X), word(X, c(2), sep = fixed("_")),as.character(X)),
         ONT_isoform = ifelse(grepl("PB",ONT_isoform),"NA",ONT_isoform))

length(unique(primer_specific$X))
length(unique(primer_specific$ONT_isoform))

annotated$ont$Norm_transcounts %>% filter(isoform %in% primer_specific$ONT_isoform) %>% 
  filter(time == 8) %>%
  mutate(Genotype = factor(ifelse(group == "CONTROL", "WT", "TG"), levels = c("WT","TG"))) %>% 
  ggplot(., aes(x = Genotype, y = log10(value))) + geom_boxplot() + geom_jitter() + 
  facet_grid(~isoform)

# overlap between exon 17 skipped and upregulated vs primer specific isoforms 
intersect(primer_specific$ONT_isoform,diff_exon17_skipped$isoform)



# 2. primer to negative control (TALONT000794296)
# also upregulated: TALONT000761829

plot_trans_exp_individual("ENSMUST00000234857.1", annotated$ont$Norm_transcounts)

annotated$ont$Norm_transcounts %>% filter(isoform == "TALONT000761829") %>%
  filter(time == 8) %>%
  ggplot(., aes(x = group, y = value)) + geom_boxplot()


mean_counts <- annotated$ont$Norm_transcounts %>% filter(associated_gene == "Bin1") %>%
  filter(time == 8) %>%
  group_by(isoform, group) %>% tally(value)

annotated$ont$Norm_transcounts %>% 
  filter(isoform == "TALONT000794296") %>% 
  mutate(Genotype = factor(ifelse(group == "CONTROL", "WT", "TG"), levels = c("WT","TG"))) %>% 
  filter(time == 8) %>%
  ggplot(., aes(x = Genotype, y = value)) + geom_boxplot() + geom_jitter() + 
  facet_grid(~isoform)  + 
  labs(x = "Genotype", y = "Expression (log10 normalised counts)")


# 3. IR event 


