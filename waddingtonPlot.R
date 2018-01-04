options(stringsAsFactors = F)
rm(list = ls())
library(ggplot2)
library(ggridges)
waddington_layout <- matrix(c(1,1,1,1,
                             1,1,2,2,
                             1,2,3,3,
                             1,2,3,4,
                             1,2,3,4,
                             1,2,3,4),
                           nrow = 6, ncol = 4, byrow = T)
waddingtonPlot(waddington_layout, do.return =T) + 
   scale_fill_gradientn(colours = c("red","orange","yellow","green","blue","cyan","purple"))                
####
waddingtonPlot <- function(layout_matrix, width_point = 1000, height_point = 60, height.scale = 5, width.ratio = 0.3, do.return = F){
  require(ggplot2)
  require(ggridges)
  # check layout_matrix is legal
  if(!all(apply(layout_matrix, 1, function(x){
    diff(x) >= 0
  }))) {
    stop("Illegal layout!")
  }

  slot.nums <- apply(layout_matrix, 1, function(x) length(unique(x)))
  slot.curves <- apply(layout_matrix, 1, function(x){
    sections <- quantile(1:width_point, probs = as.numeric(table(x)/length(x)), type = 1)
    curves <- sapply(sections, function(x) cos((1:x)*2*pi/x))
    return(unlist(curves))
  })
  
  interpolate <- quantile(1:height_point, probs = 1/(nrow(layout_matrix)-1), type = 1)
  ggData <- NULL
  height_total <- (nrow(layout_matrix)-1)*interpolate
  ratio <- 1-width.ratio
  for(i in 1:(nrow(layout_matrix)-1)){
    curve1 <- slot.curves[,i]
    curve2 <- slot.curves[,i+1]
    for(j in 0:(interpolate-1)){
      alpha <- j/interpolate
      curve <- (1-alpha)*curve1 + alpha*curve2 + 1
      height_j <- nrow(layout_matrix)*interpolate - (j + interpolate*i)
      ggData <- rbind(ggData, data.frame(x = (1:width_point)*(1-ratio*height_j/height_total) + ratio*width_point*height_j/2/height_total, 
                                         y = unname(height_j), height = curve))
    }
  }
  p <- ggplot() + geom_ridgeline(data = ggData,
                            mapping = aes(x, y, height = height + 2, group = y, fill = y),
                            scale = height.scale, color = "cyan") + theme_void()
  if(do.return){
    return(p)
  }else{
   print(p)
    return(NULL)
  }
}

