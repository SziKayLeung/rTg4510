library("stringr")
fast5Samples <- data.table::fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/sorted_nuclei/RNA/mouse/NeuN/NeuN/20240207_1339_3D_PAU69826_5fdf4c14/fastq_pass/Samples.txt", data.table = FALSE, col.names = c("Sample"))
demuxSamples <- data.table::fread("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/1_demultiplex/NeuN/processedSamples.txt", data.table = FALSE, col.names = c("Sample"))

demuxSamples <- data.table::fread("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/1_demultiplex/NeuN/processedSamples2.txt", data.table = FALSE, col.names = c("Sample"))


fast5Samples <- str_remove(fast5Samples$Sample,".fastq.gz")
setdiff(demuxSamples$Sample, fast5Samples)
setdiff(fast5Samples,demuxSamples$Sample)
missingSamples <- setdiff(fast5Samples,demuxSamples$Sample)[1:39]
write.table(missingSamples,"/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/1_demultiplex/NeuN/missingSamples.txt", quote=F, row.names=F, col.names = F)


fast5Samples <- data.table::fread("/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/sorted_nuclei/RNA/mouse/DN/20240207_1336_2C_PAU74161_a82a4529/fastq_pass/Samples.txt", data.table = FALSE, col.names = c("Sample"))
demuxSamples <- data.table::fread("/lustre/projects/Research_Project-MRC148213/lsl693/rTg4510/H_Sorted_Nuclei/1_demultiplex/DN/processedSamples.txt", data.table = FALSE, col.names = c("Sample"))
fast5Samples <- str_remove(fast5Samples$Sample,".fastq.gz")
setdiff(demuxSamples$Sample, fast5Samples)
setdiff(fast5Samples,demuxSamples$Sample)
