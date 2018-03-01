# argument
Args <- commandArgs()
nArgs <- length(Args)
species <- Args[6]
genes <- Args[7:nArgs]

# process
options(stringsAsFactors = F)

capitalize <- function(str){
  sapply(str, function(x) {
    tmp <- strsplit(x, split = "")[[1]]
    paste0(toupper(tmp[1]), paste0(tmp[-1], collapse = ""))
  })
}

geneMatch <- function(genes.1, genes.2, ignore.case = T){
  genes.2 <- unique(genes.2)
  bak.1 <- genes.1
  bak.2 <- genes.2
  
  if(ignore.case){
    genes.1 <- tolower(genes.1)
    genes.2 <- tolower(genes.2)
  }
  res <- NULL
  for(i in 1:length(genes.1)){
    matched <- bak.2[which(genes.2 == genes.1[i])]
    if(length(matched) == 0){
      res <- rbind(res, cbind(bak.1[i], bak.1[i]))
    }else{
      res <- rbind(res, 
                   cbind(bak.1[i], matched))
    }
  }
  colnames(res) <- c("Gene1", "Gene2")
  return(res)
}
alias2SymbolTable <- function(alias, species = c("Hs","Mm")){
  require(package = paste0("org.", species,".eg.db"), character.only = T)
  aliasTable <- geneMatch(alias, keys(get(paste0("org.", species,".egALIAS2EG"))))
  select(x = get(paste0("org.", species,".eg.db")), keys = aliasTable[,2], keytype = "ALIAS",columns = c("SYMBOL"))
}



if(tolower(species) %in% c("mm", "mmu","mouse","mus")) {
  species <- "Mm"
  genes <- capitalize(tolower(genes))
}else if(tolower(species) %in% c("hs","hsa","homo","human")){
  species <- "Hs"
  genes <- toupper(genes)
}

#print(genes)
res <- suppressWarnings(suppressMessages(alias2SymbolTable(alias = genes, species = species)))
write.csv(x = res, file = "gene.csv", row.names = F)
