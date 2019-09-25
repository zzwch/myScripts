'%strin%' <- function(x, table, ignore.case = T){
  unlist(sapply(x, function(i) {any(grepl(pattern = i, table, ignore.case = ignore.case))}))
}

corner <- function(mat, n = c(5,5), pos = 'top-left'){
  if(length(n) == 1) n <- rep(n, 2)
  nr <- min(n, nrow(mat))
  nc <- min(n, ncol(mat))
  ava_pos <- levels(interaction(c('top', 'middle', 'bottom'),
                                c('left', 'right', 'center'), 
                                sep = '-', lex.order = T))
  if(pos %in% ava_pos){
    if('top' %strin% pos){
      nrs <- 1
      nre <- nr+nrs-1
    }
    if('middle' %strin% pos){
      nrs <- round((nrow(mat) - nr)/2) +1
      nre <- nr+nrs-1
    }
    if('bottom' %strin% pos){
      nre <- nrow(mat)
      nrs <- nre-nr+1
    }
    
    if('left' %strin% pos){
      ncs <- 1
      nce <- nc+ncs-1
    }
    if('center' %strin% pos){
      ncs <- round((ncol(mat) - nc)/2) +1
      nce <- nc+ncs-1
    }
    if('right' %strin% pos){
      nce <- ncol(mat)
      ncs <- nce-nc+1
    }
  }else{
    stop(paste0('available position is ', paste0(ava_pos, collapse = ", ")))
  }
  
  return(mat[nrs:nre, ncs:nce])
}
