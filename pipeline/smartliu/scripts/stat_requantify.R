Args <- commandArgs()
nArgs <- length(Args)
if(nArgs != 10){
  cat("Example: Rscript stat_requantify.R summary_dir quantify_dir samples regtf_path regtf ")
  stop("please set sufficient Arguements")
}

###
# user settings
###
summary_dir <- Args[6]
quantify_dir <- Args[7]
samples <- strsplit(Args[8], split = ",")[[1]]
regtf_path <- Args[9]
regtf <- Args[10]

options(stringsAsFactors = F)

#get gene_list from gtf
genes <- sort(unique(rtracklayer::import(regtf_path)$gene_id))
genemat <- as.data.frame(matrix(0, length(genes), 0, dimnames = list(genes,NULL)))
readmat <- as.data.frame(matrix(0, length(genes), 0, dimnames = list(genes,NULL)))
for(s in samples){
  s_quantify <- file.path(quantify_dir, paste0(s, ".unmapped.",regtf,".requantify.txt"))
  tmp <- readLines(s_quantify)
  if(length(tmp) == 0) next()
  m1 <- matrix(unlist(strsplit(gsub(pattern = ":[0-9]+", "", tmp), split = "\t")), nrow = length(tmp), byrow = T)
  rownames(m1) <- m1[,1];  m1 <- m1[,-1,drop = F]; mode(m1) <- "numeric"
  m2 <- matrix(unlist(strsplit(gsub(pattern = "[0-9]+:", "", tmp), split = "\t")), nrow = length(tmp), byrow = T)
  rownames(m2) <- m2[,1];  m2 <- m2[,-1,drop = F]; mode(m2) <- "numeric"
  
  colnames(m1) <- paste0(s,"_sc", 1:ncol(m1))
  for(i in colnames(m1)) genemat[[i]] <- 0
  genemat[rownames(m1), colnames(m1)] <- m1
  
  colnames(m2) <- paste0(s,"_sc", 1:ncol(m2))
  for(i in colnames(m2)) readmat[[i]] <- 0
  readmat[rownames(m2), colnames(m2)] <- m2
}
s_regtf <- strsplit(s_quantify, "\\.")[[1]]
write.table(genemat, file = file.path(summary_dir, paste0("stat_",basename(quantify_dir),".recount_with_",s_regtf[length(s_regtf)-2],".umi.txt")), sep = "\t")
write.table(readmat, file = file.path(summary_dir, paste0("stat_",basename(quantify_dir),".recount_with_",s_regtf[length(s_regtf)-2],".read.txt")), sep = "\t")