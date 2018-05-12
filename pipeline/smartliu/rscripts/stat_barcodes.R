Args <- commandArgs()
nArgs <- length(Args)

if(nArgs != 9){
  cat("Example: Rscript stat_barcodes.R summary_dir clean_dir samples barcodes")
  stop("please set sufficient Arguements")
}

###
# user settings
###
summary_dir <- Args[6]
clean_dir <- Args[7]
samples <- strsplit(Args[8], split = ",")[[1]]
barcodes <- strsplit(Args[9], split = ",")[[1]]

options(stringsAsFactors = F)
stat_barcodes1 = NULL
stat_barcodes2 = NULL
library(rjson)
for(s in samples){
  json_file <- file.path(clean_dir, paste0(s,".barcodes.json"))
  json_data <- fromJSON(file = json_file)
  tmp <- as.data.frame(t(sapply(json_data[barcodes], unlist))[barcodes,])
  mismatch = 1:(ncol(tmp)/2)-1
  colnames(tmp) <- c(paste0("split_mis", mismatch), paste0("valid_mis", mismatch))
  rownames(tmp) <- paste0(s, "_sc", 1:length(barcodes))
  stat_barcodes1 <- rbind(stat_barcodes1, tmp)
  stat_barcodes2 <- rbind(stat_barcodes2, 
                          data.frame(split = c(rowSums(tmp[,1:(ncol(tmp)/2)]), 
                                               setNames(c(json_data$ambiguous[1], json_data$unmatched[1]),
                                                        paste0(s,c("_ambiguous","_unmatched")))
                                               ), 
                                     valid = c(rowSums(tmp[,(ncol(tmp)/2+1):ncol(tmp)]),
                                               setNames(c(json_data$ambiguous[2], json_data$unmatched[2]),
                                                        paste0(s,c("_ambiguous","_unmatched")))
                                               )
                                     )
                          )
}
write.table(stat_barcodes1, file = file.path(summary_dir, "stat_barcodes.mismatch.xls"), sep = "\t")
write.table(stat_barcodes2, file = file.path(summary_dir, "stat_barcodes.xls"), sep = "\t")
