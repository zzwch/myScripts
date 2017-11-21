############################################################################
# Produce phenotype (.cls) and dataset (.txt) files for GSEA-java analysis
############################################################################
# expr - gene expression matrix with genes in rownames, samples in colnames
#        should be data.frame or matrix.
#        only samples in rownames(annot) will be used.
# annot - sample annotation matrix with samples in rownames, annotations in colnames
#         should be data.frame or matrix
# filenames.prefix - prefix used for saved files        
# add.sys.time - add string derived from Sys.time() to the end of filenames of saved files 
write.gsea <- function(expr, annot, filename.prefix = "gsea", add.sys.time = T){
  suffix <- ifelse(add.sys.time, Sys.time(), "")
  filename.dataset <- paste(filename.prefix, suffix, "txt",sep = ".")
  expr <- expr[,rownames(annot)]
  dataset <- cbind(rownames(expr), NA, expr)
  colnames(dataset)[1:2] <- c("Name", "Description")
  write.table(dataset, file = filename.dataset, quote = F, sep = "\t", row.names = F)
  sapply(colnames(annot), function(x){
    label <- unique(annot[,x])
    filename.pheno <- paste(filename.prefix, x, suffix, "cls",sep = ".")
    write(paste(length(annot[,x]), length(label), 1, sep = " "), file = filename.pheno)
    write(paste("#", paste(label, collapse = " "), sep = " "), file = filename.pheno, append = T)
    write(paste(annot[,x], collapse = " "), file = filename.pheno,append = T)
    #unlink(filename.pheno)
    return()
  })
  return()
}

tmp <- "EP7"
write.gsea(lofExpr, subset(annot, Cluster == tmp)[,"Group",drop = F], filename.prefix = tmp)

tmp <- "EP7"
write.gsea(10*(2^lofExpr-1), subset(annot, Cluster == tmp)[,"Group",drop = F], filename.prefix = paste0("TPM.",tmp))


write.gsea(lofExpr, subset(annot, Cluster %in% c("EP6","EP7") & Group == "ht")[,"Cluster",drop = F], filename.prefix = "ht.EP6vsEP7")
