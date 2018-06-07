Args <- commandArgs()
nArgs <- length(Args)
if(nArgs != 9){
  cat("Example: Rscript stat_quantify.R summary_dir quantify_dir samples gtf_path")
  stop("please set sufficient Arguements")
}

###
# user settings
###
summary_dir <- Args[6]
quantify_dir <- Args[7]
samples <- strsplit(Args[8], split = ",")[[1]]
gtf_path <- Args[9]

options(stringsAsFactors = F)

#get gene_list from gtf
genes <- sort(unique(rtracklayer::import(gtf_path)$gene_id))
umimat <- as.data.frame(matrix(0, length(genes), 0, dimnames = list(genes,NULL)))
readmat <- as.data.frame(matrix(0, length(genes), 0, dimnames = list(genes,NULL)))
mapmat <- NULL
for(s in samples){
  s_quantify <- file.path(quantify_dir, paste0(s, ".quantify.txt"))
  tmp <- readLines(s_quantify)
  if(length(tmp) == 0) next()
  m1 <- matrix(unlist(strsplit(gsub(pattern = ":[0-9]+", "", tmp), split = "\t")), nrow = length(tmp), byrow = T)
  rownames(m1) <- m1[,1];  m1 <- m1[,-1]; mode(m1) <- "numeric"
  m2 <- matrix(unlist(strsplit(gsub(pattern = "[0-9]+:", "", tmp), split = "\t")), nrow = length(tmp), byrow = T)
  rownames(m2) <- m2[,1];  m2 <- m2[,-1]; mode(m2) <- "numeric"
  
  colnames(m1) <- paste0(s,"_sc", 1:ncol(m1))
  for(i in colnames(m1)) umimat[[i]] <- 0
  umimat[rownames(m1), colnames(m1)] <- m1
  
  colnames(m2) <- paste0(s,"_sc", 1:ncol(m2))
  for(i in colnames(m2)) readmat[[i]] <- 0
  readmat[rownames(m2), colnames(m2)] <- m2
  
  m3 <- read.table(file = file.path(quantify_dir, paste0(s, ".quantify.sam.flagstat")), header = F, sep = "\t", row.names = 1)
  colnames(m3) <- paste0(s,"_sc", 1:ncol(m3))
  if(is.null(mapmat)) mapmat <- m3
  else{
    if(all(sort(rownames(mapmat)) == sort(rownames(m3)))){
      mapmat <- cbind(mapmat, m3[rownames(mapmat),])
    }else{
      m3_rn <- sort(unique(c(rownames(mapmat), rownames(m3))))
      tmp <- as.data.frame(matrix(0, length(m3_rn), 0, dimnames = list(m3_rn,NULL)))
      for(i in m3_rn) mapmat[[i]] <- 0
      mapmat[rownames(m3), colnames(m3)] <- m3
      mapmat[rownames(mapmat), colnames(mapmat)] <- mapmat
    }
  }
}
umimat <- cbind(gene_symbol = rownames(umimat), umimat)
readmat <- cbind(gene_symbol = rownames(readmat), readmat)
mapmat <- cbind(map_flag = rownames(mapmat), mapmat)
write.table(umimat, file = file.path(summary_dir, paste0("stat_",basename(dirname(quantify_dir)),".", basename(quantify_dir),".umi.txt")), sep = "\t", row.names = F)
write.table(readmat, file = file.path(summary_dir, paste0("stat_",basename(dirname(quantify_dir)),".", basename(quantify_dir),".read.txt")), sep = "\t", row.names = F)
write.table(mapmat, file = file.path(summary_dir, paste0("stat_",basename(dirname(quantify_dir)),".", basename(quantify_dir),".flag.txt")), sep = "\t", row.names = F)