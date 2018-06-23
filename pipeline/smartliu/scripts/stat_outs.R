Args <- commandArgs()
nArgs <- length(Args)
if(nArgs != 7){
  cat("Example: Rscript stat_outs.R summary_dir outs_dir")
  stop("please set sufficient Arguements")
}

###
# user settings
###
summary_dir <- Args[6]
outs_dir <- Args[7]
# summary_dir <- "../summary/"
# outs_dir <- "./"

options(stringsAsFactors = F)

# pool mix estimate
require(mixtools)
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma) * length(x)
}


require(ggplot2)
require(reshape2)
require(ggthemr)
require(gridExtra)
ggthemr(palette = "flat", layout = "clear", spacing = 1, type = "outer")
#################################
# get barcodes splitting data
##############
# by pool
barcodes_split <- as.data.frame(read.delim(file.path(summary_dir,"stat_barcodes.txt"), header = T, check.names = F))
pool_n <- nrow(barcodes_split)/97
barcodes_split$Sample <- gsub(pattern = "_sc[0-9]+|_unmatched$", "", barcodes_split$cell_id)
barcodes_split$Match <- ifelse(grepl("_sc[0-9]+$",barcodes_split$cell_id), "matched", "unmatched")

stat_split <- aggregate(barcodes_split[,c("split","valid")], list(Sample = barcodes_split$Sample, Match = barcodes_split$Match), sum)
stat_split <- rbind(stat_split, cbind(aggregate(stat_split[,c("split","valid")], list(Sample = stat_split$Sample), sum), Match = "whole"))
stat_split$valid.ratio <- stat_split$valid/stat_split$split
stat_split$match.ratio <- apply(stat_split, 1, function(x){
  tmp <- stat_split$Sample == x["Sample"]
  stat_split[tmp & stat_split$Match == "matched", "split"]/stat_split[tmp & stat_split$Match == "whole", "split"]
})

stat_split$pool.magnitudes <- apply(stat_split, 1, function(x){
  tmp <- barcodes_split[barcodes_split$Sample == x["Sample"] & barcodes_split$Match == "matched", "split"]
  mixmdl <- normalmixEM(log10(tmp+1), k = 2)
  ind_valid <- which(mixmdl$posterior[,which.max(mixmdl$mu)] > 0.6)
  mean(log10(tmp[ind_valid]+1))
})
stat_split$pool.variance <- apply(stat_split, 1, function(x){
  tmp <- barcodes_split[barcodes_split$Sample == x["Sample"] & barcodes_split$Match == "matched", "split"]
  mixmdl <- normalmixEM(log10(tmp+1), k = 2)
  ind_valid <- which(mixmdl$posterior[,which.max(mixmdl$mu)] > 0.6)
  var(log10(tmp[ind_valid]+1))
})

# plot barcode counts
ggData <- melt(stat_split[stat_split$Match != "whole", c("Sample","Match","split","valid")],
               id.vars = c("Sample","Match"), variable.name = "Valid",value.name = c("Count"))
p_pool_count <- ggplot(ggData) + 
  geom_bar(aes(Valid, Count,fill = Match), position = position_stack(), stat = "identity") + 
  facet_wrap(~Sample, nrow = 1) + theme(legend.position = "top") +
  labs(title = "Counting barcode-match reads after 96-barcodes splitting", 
       subtitle = "by Pool with consideration of barcode mismatch and read length validation")
# plot barcode ratio
p_pool_ratio <- ggplot(stat_split, mapping = aes(Sample, valid.ratio, color = Match)) + 
  geom_point() + geom_line() +
  theme(legend.position = "top") + ylim(c(0,1)) +
  labs(title = "Ratio of length-valid reads in barcode-match reads", 
       subtitle = "by Pool with consideration of barcode mismatch")
# plot pool magnitudes
p_pool_magnitudes <- ggplot(stat_split, mapping = aes(Sample, pool.magnitudes)) +
  geom_point() + geom_line() + 
  theme(legend.position = "top") +
  labs(title = "Magnitude of barcode-match reads (in log10 scale)", 
       subtitle = "by Pool")
# plot pool variance
p_pool_variance <- ggplot(stat_split, mapping = aes(Sample, pool.variance)) +
  geom_point() + geom_line() +
  theme(legend.position = "top") +
  labs(title = "Variance of barcode-match reads in pooled cells", 
       subtitle = "by Pool")
# save outputs
pdf(file = file.path(outs_dir, "outs_pool.pdf"), width = 14+ 1*pool_n, height = 10)
grid.arrange(p_pool_count, p_pool_ratio, p_pool_magnitudes, p_pool_variance, nrow = 2)
dev.off()
write.table(stat_split, file = file.path(outs_dir, "outs_pool.xls"), sep = "\t", row.names = F)

##############
# by cell
barcodes_mismatch <- as.data.frame(read.delim(file.path(summary_dir,"stat_barcodes.mismatch.txt"), header = T, check.names = F))
stat_barcodes <- barcodes_mismatch
stat_barcodes$Sample <- gsub(pattern = "_sc[0-9]+$", "", stat_barcodes$cell_id)
stat_barcodes$Barcode <- as.numeric(gsub(pattern = ".*_sc([0-9]+)$", "\\1", stat_barcodes$cell_id))

stat_barcodes$split <- rowSums(barcodes_mismatch[,grepl("split_mis", colnames(barcodes_mismatch))])
stat_barcodes$valid <- rowSums(barcodes_mismatch[,grepl("valid_mis", colnames(barcodes_mismatch))])
stat_barcodes$valid.ratio <- stat_barcodes$valid/stat_barcodes$split
for(i in (1:(ncol(barcodes_mismatch)/2))-1){
  stat_barcodes[[paste0("valid.ratio.mis", i)]] <- barcodes_mismatch[[paste0("valid_mis", i)]]/barcodes_mismatch[[paste0("split_mis", i)]]
}
# plot counts of split and valid reads
ggData <- melt(stat_barcodes[,c(1,grep(pattern = "_mis", colnames(stat_barcodes)))], id.vars = c("cell_id"), variable.name = "type", value.name = "Count")
ggData$Sample <- gsub(pattern = "_sc[0-9]+", "", ggData$cell_id)
ggData$Barcode <- as.numeric(gsub(pattern = ".*_sc([0-9]+)", "\\1", ggData$cell_id))
ggData$Group <- gsub(pattern = "_mis[0-9]", "", ggData$type)
ggData$mismatch <- gsub(pattern = ".*_mis", "", ggData$type)

p_barcode_count <- ggplot(ggData) +  
  geom_bar(mapping = aes(Barcode, Count, fill = mismatch, group = Group), 
           position = "dodge", stat = "identity", width = 0.6) +
  facet_wrap(~Sample, ncol = 1) + theme(legend.position = "top") +
  scale_x_continuous(minor_breaks = seq(0,96, 5)) +
  theme(panel.grid.minor.x = element_line(colour = "grey", linetype = "dashed")) +
  labs(title = "Counting barcode-match reads after 96-barcodes splitting", 
       subtitle = "by each cell with consideration of barcode mismatch and read length validation")

# plot ratio of valid versus split reads
ggData <- melt(stat_barcodes[,c(1,grep(pattern = "valid.ratio", colnames(stat_barcodes)))], id.vars = c("cell_id"), variable.name = "type", value.name = "Ratio")
ggData$Sample <- gsub(pattern = "_sc[0-9]+", "", ggData$cell_id)
ggData$Barcode <- as.numeric(gsub(pattern = ".*_sc([0-9]+)", "\\1", ggData$cell_id))
#ggData[is.na(ggData)] <- 0
p_barcode_ratio <- ggplot(ggData, mapping = aes(Barcode, Ratio, color = type)) +  
  geom_line(size = 1) + geom_point(size = 1) +
  facet_wrap(~Sample, ncol = 1) +  theme(legend.position = "top") +
  scale_x_continuous(minor_breaks = seq(0,96, 5)) +
  theme(panel.grid.minor.x = element_line(colour = "grey", linetype = "dashed")) +
  labs(title = "Ratio of length-valid reads in barcode-match reads", 
       subtitle = "by each cell with consideration of barcode mismatch")

# save outputs
pdf(file = file.path(outs_dir, "outs_barcodes.pdf"), width = 15, height = 5* pool_n)
grid.arrange(p_barcode_count, p_barcode_ratio, ncol = 2)
dev.off()
write.table(stat_barcodes, file = file.path(outs_dir, "outs_barcodes.xls"), sep = "\t", row.names = F)

#################################
# get mapping data
m_files <- list.files(path = summary_dir, pattern = "^stat_mapping_to_")
m_refs <- gsub(pattern = "^stat_mapping_to_|\\.(re)?(count|mapping).*|\\.txt","",m_files, perl = T)
for(i in sort(unique(m_refs))){
  ref_stat <- NULL
  ref_mapping_file <- paste0("stat_mapping_to_", i, ".txt")
  
  ref_quantify_files <- setdiff(m_files[m_refs == i], ref_mapping_file)
  ref_gtfs <- gsub(".*\\.count_with_|\\.(flag|umi|read)\\.txt", replacement = "", ref_quantify_files[grep("\\.count_with_",ref_quantify_files)])
  ref_regtfs <- gsub(".*\\.recount_with_|\\.(flag|umi|read)\\.txt", replacement = "", ref_quantify_files[grep("\\.recount_with_",ref_quantify_files)])
  for(j in sort(unique(ref_gtfs))){
    ref_gtf_flag <- paste0("stat_mapping_to_", i, ".count_with_",j,".flag.txt")
    ref_gtf_read <- paste0("stat_mapping_to_", i, ".count_with_",j,".read.txt")
    ref_gtf_umi <- paste0("stat_mapping_to_", i, ".count_with_",j,".umi.txt")
    # save outputs
    if(!file.copy(from = file.path(summary_dir, ref_gtf_umi), to = file.path(outs_dir, paste0(i,".",j,".umi.txt")), overwrite = T)) stop("file.copy error")
    
    flag <- as.data.frame(t(read.delim(file.path(summary_dir, ref_gtf_flag), header = T, check.names = F, row.names = 1)))
    colnames(flag) <- plyr::mapvalues(x = colnames(flag), paste0("flag_", c(0,4,16,256,272)), c("mapped_forward", "unmapped", "mapped_reverse", "mapped_notprimary_forward", "mapped_notprimary_reverse"))
    read <- as.data.frame(read.delim(file.path(summary_dir, ref_gtf_read), header = T, check.names = F, row.names = 1))
    umi <- as.data.frame(read.delim(file.path(summary_dir, ref_gtf_umi), header = T, check.names = F, row.names = 1))
    ref_gtf_stat <- cbind(flag, gene_reads = colSums(read), gene_umis = colSums(umi), 
                          spike_reads = colSums(read[grep("^ERCC-|^RGC-", rownames(read)),,drop = F]),
                          spike_umis = colSums(umi[grep("^ERCC-|^RGC-", rownames(umi)),,drop = F]))
    ref_gtf_stat$barcode <- as.numeric(gsub(".*_sc", "", rownames(ref_gtf_stat)))
    ref_gtf_stat$sample <- gsub("_sc[0-9]+$", "", rownames(ref_gtf_stat))
    ref_gtf_stat$mapping_rate <- rowSums(ref_gtf_stat[,c("mapped_forward","mapped_reverse")])/rowSums(ref_gtf_stat[,c("mapped_forward","mapped_reverse","unmapped")])
    ref_gtf_stat$umi_rate <- ref_gtf_stat$gene_umis/ref_gtf_stat$gene_reads
    ref_gtf_stat$spike_ratio_read <- ref_gtf_stat$spike_reads/ref_gtf_stat$gene_reads
    ref_gtf_stat$spike_ratio_umi <- ref_gtf_stat$spike_umis/ref_gtf_stat$gene_umis
    
    ##############
    # plot count
    ggData <- melt(ref_gtf_stat[,c("sample", "barcode", "gene_reads", "gene_umis")], 
                   id.vars = c("sample","barcode"), variable.name = "type", value.name = "count")
    p_cell_count <- ggplot(ggData) +  
      geom_bar(mapping = aes(barcode, count, fill = type, group = type), 
               position = "dodge", stat = "identity", width = 0.6) +
      facet_wrap(~sample, ncol = 1) + theme(legend.position = "top") +
      scale_x_continuous(minor_breaks = seq(0,96, 5)) +
      theme(panel.grid.minor.x = element_line(colour = "grey", linetype = "dashed")) +
      labs(title = "Counting reads/UMIs mapped to gene exons", 
           subtitle = "by each cell")
    
    # plot ratio
    ggData <- melt(ref_gtf_stat[,c("sample", "barcode", "mapping_rate", "umi_rate", "spike_ratio_read", "spike_ratio_umi")], 
                   id.vars = c("sample","barcode"), variable.name = "type", value.name = "ratio")
    p_cell_ratio <- ggplot(ggData, mapping = aes(barcode, ratio, color = type)) +  
      geom_line(size = 1) + geom_point(size = 1) +
      facet_wrap(~sample, ncol = 1) +  theme(legend.position = "top") +
      scale_x_continuous(minor_breaks = seq(0,96, 5)) +
      theme(panel.grid.minor.x = element_line(colour = "grey", linetype = "dashed")) +
      labs(title = "Ratio of mapped reads, unique UMIs, spkie_reads and spike_UMIs", 
           subtitle = "by each cell")
    
    pdf(file = file.path(outs_dir, paste0("outs_",i,".",j,".pdf")), width = 15, height = 5* pool_n)
    grid.arrange(p_cell_count, p_cell_ratio, ncol = 2)
    dev.off()
    ###########
    colnames(ref_gtf_stat) <- paste0(i,"|",j,"|", colnames(ref_gtf_stat))
    
    samples <- c(rownames(ref_stat), rownames(ref_gtf_stat))
    annots <- c(colnames(ref_stat), colnames(ref_gtf_stat))
    ref_stat <- as.data.frame(matrix(NA, length(samples), length(annots), dimnames = list(samples, annots)))
    ref_stat[rownames(ref_stat), colnames(ref_stat)] <- ref_stat
    ref_stat[rownames(ref_gtf_stat), colnames(ref_gtf_stat)] <- ref_gtf_stat
    
    
  }
  ref_mapping <- as.data.frame(read.delim(file.path(summary_dir,ref_mapping_file), header = T, check.names = F))
  colnames(ref_mapping) <- paste0(i,"|",tolower(colnames(ref_mapping)))
  ref_stat <- cbind(ref_stat, ref_mapping)
  for(m in 1:nrow(ref_mapping)){
    ref_stat[paste0(ref_mapping[m,1], "_sc", 1:96),colnames(ref_mapping)] <- ref_mapping[m,]
  }
  
  for(k in sort(unique(ref_regtfs))){
    ref_regtf_map <- paste0("stat_mapping_to_", i, ".remapping_to_",k,".txt")
    ref_regtf_read <- paste0("stat_mapping_to_", i, ".recount_with_",k,".read.txt")
    ref_regtf_umi <- paste0("stat_mapping_to_", i, ".recount_with_",k,".umi.txt")
    # save outputs
    file.copy(from = file.path(summary_dir, ref_regtf_umi), to = file.path(outs_dir, paste0(i,".requantify.",k,".umi.txt")), overwrite = T)
    
  }
  # save outputs
  
  write.table(cbind(cell_id =rownames(ref_stat), ref_stat), file = file.path(outs_dir, paste0("outs_",i,".xls")), sep = "\t", row.names = F)
  
}
