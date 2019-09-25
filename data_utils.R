'%strin%' <- function(x, y){
  
}

corner <- function(x, n = c(5,5), pos = 'top-left') {
  if(length(n) = 1) n <- rep(n, 2)
  nr <- min(n, nrow(x))
  nc <- min(n, ncol(x))
  if(pos %in% levels(interaction(c('top', 'middle', 'bottom'),
                   c('left', 'right', 'center'), 
                   sep = '-', lex.order = T))){
      
  }

  x[1:min(n, nrow(x)), 1:min(n, ncol(x)]
}
