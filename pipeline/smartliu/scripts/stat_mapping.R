Args <- commandArgs()
nArgs <- length(Args)
if(nArgs != 8){
  cat("Example: Rscript stat_mapping.R summary_dir mapping_dir samples")
  stop("please set sufficient Arguements")
}

###
# user settings
###
summary_dir <- Args[6]
mapping_dir <- Args[7]
samples <- strsplit(Args[8], split = ",")[[1]]

options(stringsAsFactors = F)

#Hisat2_log: mapping summary
alignStat <- matrix(NA,length(samples),6)
rownames(alignStat) <- samples
colnames(alignStat) <- c("Sample","Input_reads", "Mapped_reads", "Mapping_rate", "Unique_mapped_reads", "Unique_mapping_rate")
for(s in samples){
  align_log <- file.path(mapping_dir, paste0(s,".align_summary.txt"))
  tmp <- try(read.table(file = align_log, skip = 1, nrows = 4, sep = "\t"), silent = F)
  if(class(tmp)=="try-error") next
  tmp <- as.numeric(apply(tmp, 1, function(x) gsub(" |\\(.*","",x)))
  alignStat[s,] <- c(s, tmp[1], tmp[1]-tmp[2], (tmp[1] - tmp[2])/tmp[1], tmp[3], tmp[3]/tmp[1])
}

write.table(alignStat, file = file.path(summary_dir, paste0("stat_", basename(mapping_dir),".txt")), sep = "\t", row.names = F)
