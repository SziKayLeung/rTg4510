class.files$targ_qc %>% filter(isoform %in% class.files$targ_filtered$isoform)

junc.files.targ_qc <- data.table::fread(paste0(dirnames$targ_root, "/3_sqanti3/reRunValidation/all_iso_ont_collapsed_junctions.txt"), header = T)

class.files$targ_qc <- class.files$targ_qc %>% filter(isoform %in% class.files$targ_filtered$isoform) 
table(class.files$targ_qc$RNASeq_supported)
table(class.files$targ_qc$within_50_cage)
table(class.files$targ_qc$RNASeq_supported)
table(class.files$targ_qc$polyA_motif_found)
table(class.files$targ_qc$within_polya_site)

junc.files.targ_qc <- data.frame(junc.files.targ_qc %>% filter(isoform %in% class.files$targ_filtered$isoform))

nrow(junc.files.targ_qc[junc.files.targ_qc$total_coverage_unique >= 1 ,])
nrow(junc.files.targ_qc)
table(junc.files.targ_qc[junc.files.targ_qc$total_coverage_unique == 0,"junction_category"])
table(junc.files.targ_qc[junc.files.targ_qc$total_coverage_unique == 0,"splice_site"])
table(junc.files.targ_qc[junc.files.targ_qc$total_coverage_unique == 0,"canonical"])
length(unique(junc.files.targ_qc[junc.files.targ_qc$total_coverage_unique == 0,"isoform"]))

length(unique(junc.files.targ_qc$isoform))


CoveragePerIsoform <- junc.files.targ_qc %>% 
  dplyr::group_by(as.factor(isoform)) %>% 
  dplyr::summarise(sum_total_coverage_unique = sum(total_coverage_unique), median_total_coverage_unique = median(total_coverage_unique))

median(CoveragePerIsoform$sum_total_coverage_unique)
median(CoveragePerIsoform$median_total_coverage_unique)
