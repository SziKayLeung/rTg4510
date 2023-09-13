
### ------------- Table 1 -------------------
mTable1 <- Merged_gene_class_df[c("Total.Number.of.Transcripts","Total.Number.of.Known.Transcripts","Total.Number.of.Novel.Transcripts","Total.Number.of.Coding.Transcripts",
                       "Total.Number.of.Transcripts.with.canonical.splice.sites","Total.Number.of.Transcripts.with.non.canonical.splice.sites",
                       "A5A3","ES","IR",
                       "Number.of.Transcripts.with.A5A3.Events","Number.of.Transcripts.with.ES.Events","Number.of.Transcripts.with.IR.Events"),]
Table1 <- data.table::transpose(mTable1)
colnames(Table1) <- c("All","Known","Novel","Coding","Canonical","Non canonical","A5A3","ES","IR","A5A3T","EST","IRT")
rownames(Table1) <- colnames(mTable1)
Table1$NonCoding <-(Table1$All-Table1$Coding)
Table1 <- Table1 %>% mutate(A5A3F = paste0(A5A3T," (",A5A3,")"), 
                  ESF = paste0(EST," (",ES,")"),
                  IRF = paste0(IRT," (",IR,")")) 
Table1 <- Table1 %>% select(All,Known,Novel,Coding,NonCoding,A5A3F,ESF,IRF)


### ---------- Table 2 --------------------
class.files$targ_filtered[class.files$targ_filtered$isoform %in% TargetedDESeq$ontResTranAnno$wald$anno_res$isoform,c("isoform","structural_category","subcategory")]
