homoconv <- function(msigdb.dir = NULL, out.dir = NULL, use.biomaRt = T, hcop.file = NULL, try.most = T){
  if(!dir.exists(out.dir)){
    stop("out.dir is invalid!")
  }
  options(stringsAsFactors = F)
  if(is.null(msigdb.dir)){
    stop("\nPlease set MsigDB dir which includes genesets files")
  }

  
  msigdb.files <- list.files(path = msigdb.dir)
  msigdb.entrez <- grep(pattern = "entrez.gmt", x = msigdb.files, ignore.case = T, value = T)
  msigdb.symbols <- grep(pattern = "symbols.gmt", x = msigdb.files, ignore.case = T, value = T)
  
  # entrez
  genesets.entrez <- lapply(msigdb.entrez,function(x){
    read.gmt(gmtfile = paste0(msigdb.dir,"/",x))
  })
  names(genesets.entrez) <- msigdb.entrez
  all.entrez <- Reduce(union, lapply(genesets.entrez, function(x) unlist(x$gmt)))
  # symbol
  genesets.symbols <- lapply(msigdb.symbols,function(x){
    read.gmt(gmtfile = paste0(msigdb.dir,"/",x))
  })
  names(genesets.symbols) <- msigdb.symbols
  all.symbols <- Reduce(union, lapply(genesets.symbols, function(x) unlist(x$gmt)))
  
  if(use.biomaRt){
    if(is.null(hcop.file)){
      require(biomaRt)
      human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      all.symbols.mapping <- getLDS(attributes = c("hgnc_symbol"),
                                    filters = "hgnc_symbol", values = all.symbols, mart = human,
                                    attributesL = c("mgi_symbol"), martL = mouse)
      colnames(all.symbols.mapping) <- c("human","mouse")
      
      all.entrez.mapping <- getLDS(attributes = c("entrezgene"),
                                   filters = "entrezgene", values = all.entrez, mart = human,
                                   attributesL = c("entrezgene"), martL = mouse)
      colnames(all.entrez.mapping) <- c("human","mouse")
      all.entrez.mapping[["human"]] <- as.character(all.entrez.mapping[["human"]])# may be speed up when %in%
      all.entrez.mapping[["mouse"]] <- as.character(all.entrez.mapping[["mouse"]])#
      
    }else{
      stop("Please set hcop.file = NULL when use.biomaRt = TRUE, otherwise confusing!")
    }
  }else{
    if(is.null(hcop.file)){
      stop("Please set hcop.file to the homolog file when use.biomaRt = FALSE, otherwise confusing!\nPlease download mapping file from HGNC Comparison of Orthology Predictions (HCOP). http://www.genenames.org/cgi-bin/hcop")
    }else{
      stop("Sorry, havn't been done! contact lizc07@vip.qq.com")
      
      ## to be continued
      tmp <- try(hcop <- read.table(file = hcop.file, header = T, sep = "\t", row.names = NULL, quote = ""))
      if(class(tmp) == "try-error"){
        stop("\nCannot read your provided file!\nPlease download mapping file from HGNC Comparison of Orthology Predictions (HCOP). http://www.genenames.org/cgi-bin/hcop")
      }else{
        all.entrez.mapping <- unique(hcop[which(hcop$human_entrez_gene %in% all.entrez), c("human_entrez_gene", "mouse_entrez_gene", "mouse_symbol")])
        all.symbols.mapping <- unique(hcop[which(hcop$human_symbol %in% all.symbols), c("human_symbol", "mouse_entrez_gene", "mouse_symbol")])
        
       
      }
      
    }
    
  }
  
  lapply(names(genesets.entrez), function(x){
    gmt <- genesets.entrez[[x]][["gmt"]]
    desc <- genesets.entrez[[x]][["desc"]]
    for(i in names(gmt)){
      gmt[[i]] <- unique(all.entrez.mapping[["mouse"]][which(all.entrez.mapping[["human"]] %in% gmt[[i]])])
    }
    write.gmt(gmt, desc, paste0(out.dir,"/",x))
  })
  lapply(names(genesets.symbols), function(x){
    gmt <- genesets.symbols[[x]][["gmt"]]
    desc <- genesets.symbols[[x]][["desc"]]
    for(i in names(gmt)){
      gmt[[i]] <- unique(all.symbols.mapping[["mouse"]][which(all.symbols.mapping[["human"]] %in% gmt[[i]])])
    }
    write.gmt(gmt, desc, paste0(out.dir,"/",x))
  })
  
  return()
}

read.gmt <- function(gmtfile = NULL){
  gmt <- list()
  desc <- list() 
  for(i in readLines(con = gmtfile)){
    tmp <- strsplit(i, split = "\t")[[1]]
    gmt[[tmp[1]]] <- tmp[-c(1,2)] 
    desc[[tmp[1]]] <- tmp[2]
  }
  return(list(gmt = gmt, desc = desc))
}

write.gmt <- function(gmt.list = NULL, desc.list = NULL, gmtfile = NULL){
  text <- NULL
  for(i in names(gmt.list)){
    text <- c(text, paste(c(i, desc.list[[i]], gmt.list[[i]]), collapse = "\t"))
  }
  writeLines(text = text, con = gmtfile)
  return()
}
