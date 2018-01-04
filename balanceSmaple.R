# x - a vetor with distinct categories.
# each.size - sampling size of each category
# seed - random seed
# exact - exactly same sampling size 
# low.cutoff - category whose number is less than low.cutoff will be removed
# replace - sample with replace or not.
balanceSample <- function(x, each.size = 50, seed = 666, exact = T, low.cutoff = 0.1, replace = T){
  x.tab <- table(x)
  x.tab <- x.tab[x.tab/each.size > low.cutoff]
  if(!replace & any(each.size > x.tab)){
    stop("Please set replace=T if each.size is too large")
  }
  ret.ind <- NULL
  set.seed(seed)
  if(exact){
    for(i in names(x.tab)){
      ret.ind <- c(ret.ind, sample(which(x == i), size = each.size, replace = replace))
    }
  }else{
    x.prob <- 1/x.tab[x]
    x.prob[is.na(x.prob)] <- 0
    ret.ind <- sample(1:length(x), length(x.tab)*each.size, replace = replace, prob = x.prob)
  }
  return(sort(ret.ind))
}
