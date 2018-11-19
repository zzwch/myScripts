#' this is a custom function
#' @param path specify a tar archive file downloaded from GEO, or a folder including files extracted from that tar achrvie
#' @param sep try comma if encounter a problem
#' @param show_progress boolean. show the progress bar or not.
#' @return data matrix
#' @export
readGSERAW <- function(path = NULL, # specify a tar archive file downloaded from GEO
                                    # or a folder including files extracted from that tar achrvie
                       sep = "\t", # try comma if encounter a problem
                       show_progress = T # show the progress bar or not.
                       )
{
  if(is.null(path)){
    stop("path must be specified!")
  }
  
  if(dir.exists(path)){
    gsmfiles <- list.files(path, full.names = T)
  }else{
    if(grepl(pattern = "\\.tar$|[.]", x = path, ignore.case = T, perl = T)){
      gsmfiles <- untar(tarfile = path, list = T)
      untar(tarfile = path)
    }else {
      stop("your path setting is unsupported!")
    }
  }
  gsmfiles <- gsmfiles[grepl(pattern = "/GSM", x = gsmfiles)] # in case of files copy from MAC OS
  if(show_progress) {
    pb.use <- ifelse(Sys.info()['sysname'] != 'Darwin' & require('progress'), T, F)
    n = length(gsmfiles) +1
    i <- 0
    if(pb.use){
      pb <- progress_bar$new(
        format = "Reading [:bar] :percent eta: :eta",
        total = n, clear = FALSE, width= 60)
      pb$tick(i)
    }else{
      pb <- txtProgressBar(min = 0, max = n, char = "#", title = "Reading files...",style = 3)
      setTxtProgressBar(pb, i)
    }
  } 
  # try loading one file to get gene names which is usually identical among data files.
  tmp <- read.table(file = gsmfiles[1], header = F, row.names = 1, sep = sep, stringsAsFactors = F)
  # initialize data matrix. Predefined dimesion should be more efficient than `cbind` manipulation
  gsmdata <- matrix(data = NA, 
                    nrow = nrow(tmp),
                    ncol = length(gsmfiles),
                    dimnames = list(rownames(tmp), basename(gsmfiles)))
  i <- i+1
  for(f in gsmfiles){
    if(show_progress){
      if(pb.use) pb$tick(i) else setTxtProgressBar(pb, i)
    } 
    tmp <- read.table(file = f, header = F, row.names = 1, sep = sep, stringsAsFactors = F)
    gsmdata[rownames(tmp), basename(f)] <- tmp[,1]
    i <- i+1
  }
  if(show_progress & !pb.use) close(pb)
  return(gsmdata)
}

