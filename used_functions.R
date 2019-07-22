# print a gene vector out as nrow * ncol
write.csv(formatGeneList(consensus_genes, nrow = 25), file = "clipboard", row.names = F, quote = F)
formatGeneList <- function(x, nrow = NULL, ncol = NULL){
  if(is.null(nrow) & is.null(ncol)) stop('specify nrow or ncol by yourself')
  if(!(is.null(nrow) | is.null(ncol))) stop('do not specify nrow and ncol concurrently')
  if(is.null(nrow)){
    nrow <- ceiling(length(x)/ncol)
  }
  if(is.null(ncol)){
    ncol <- ceiling(length(x)/nrow)
  }

  res <- list()
  for(i in 1:nrow){
    ind_start <- (i-1)*ncol +1
    ind_end <- i*ncol
    res[[as.character(i)]] <- paste0(x[ind_start:ind_end], collapse = ", ")
  }
  Reduce(f = function(x, y) {paste(x, y, sep = ",\n")}, x = res)
}

