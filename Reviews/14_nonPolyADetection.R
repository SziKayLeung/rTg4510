## ---------- Iso-Seq whole dataset -----------------

IsoSeq_polyAErrorRate <- function(limaPath, RefinePath){
  polyA <- lapply(list.files(RefinePath, pattern = "flnc.report.csv", full.names = T), function(x) read.csv(x))
  polyACounts <- lapply(polyA, function(x) nrow(x))
  names(polyACounts) <- list.files(RefinePath, pattern = "flnc.report.csv")
  
  lima <- lapply(list.files(path = limaPath, pattern = ".fl.lima.counts", full.names = T), function(x) read.table(x, header = T))
  names(lima) <- list.files(path = limaPath, pattern = ".fl.lima.counts")
  
  polyAFinal <- as.data.frame(do.call(rbind,polyACounts)) %>% tibble::rownames_to_column(., var = "file") %>% 
    mutate(sample = word(file,c(1),sep=fixed("."))) %>%
    `colnames<-`(c("file", "polyA", "sample")) %>% 
    select(sample, polyA)
  
  limaFinal <- as.data.frame(do.call(rbind,lima)) %>% tibble::rownames_to_column(., var = "file") %>% 
    mutate(sample = word(file,c(1),sep=fixed("."))) %>%
    select(sample, Counts) %>%
    `colnames<-`(c("sample", "lima"))
  
  PolyAErrorRate = merge(polyAFinal, limaFinal, by = "sample") %>%
    mutate(nonPolyA = lima - polyA, percNonPolyA = nonPolyA/lima * 100) 
  
  return(PolyAErrorRate)
}


IsoSeq_polyAErrorRate(limaPath = paste0(dirnames$glob_root,"/1_isoseq3/2_lima/2_lima"), 
                      RefinePath = paste0(dirnames$glob_root,"/1_isoseq3/3_refine"))

IsoSeqTargetedPolyA <- IsoSeq_polyAErrorRate(
  limaPath = paste0(dirnames$targ_iso_root,"/2_lima"), 
  RefinePath = paste0(dirnames$glob_root,"/1_isoseq3/3_refine")
)


mean(PolyAErrorRate$percNonPolyA)
mean(PolyAErrorRate$nonPolyA)
mean(PolyAErrorRate$lima)

###
ONTTargetedPolyA <- rbind(
  read.table(paste0(dirnames$targ_ont_root,"/2_demultiplex/Batch2_Demultiplex/All_cutadapt.log"), sep = ":") %>% mutate(Batch = 2),
  read.table(paste0(dirnames$targ_ont_root,"/2_demultiplex/Batch3_Demultiplex/All_cutadapt.log"), sep = ":")  %>% mutate(Batch = 3)
)

ONTTargetedPolyA <- ONTTargetedPolyA  %>% filter(V2 == "Reads with adapters") %>% 
  mutate(perc = as.numeric(word(word(V3,c(2),sep=fixed("(")),c(1),sep=fixed("%"))),
         kept = as.numeric(gsub(",", "", word(V3,c(1),sep=fixed(" (")))),
         sample = word(V1,c(1),sep=fixed("_"))) %>% 
  filter(sample != "BC10") %>% 
  mutate(totalReads = as.integer(kept / perc * 100)) %>% 
  group_by(sample, Batch) %>% 
  dplyr::summarise(n_totalReads = sum(totalReads), n_kept = sum(kept)) %>% 
  as.data.frame() %>%
  mutate(n_removed = n_totalReads - n_kept,
         perc_removed = n_removed/n_totalReads * 100)

mean(ONTTargetedPolyA$n_totalReads)
mean(ONTTargetedPolyA$n_removed)
mean(ONTTargetedPolyA$perc_removed)

