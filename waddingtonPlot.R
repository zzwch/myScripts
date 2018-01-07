options(stringsAsFactors = F)
rm(list = ls())
library(ggplot2)
library(ggridges)
waddington_layout <- list(row = c(1,1,2,3,4,4), 
                          row.scale = c(1, 1, 1, 1, 3),
                          col = list(c(1,1),
                                     c(1,1),
                                     c(1,1, 0.8,1.2),
                                     c(1.2,1, 0.9,0.9, 1.7,2.3),
                                     c(1.2,1, 0.9,0.9, 1.1,0.9, 1,1),
                                     c(1,1, 1,1, 1,1, 1,1)))
waddingtonPlot(waddington_layout, line.color = "white", line.size = 0.1, line.alpha = 0.1, line.type = 1,
               width = 1000, wave.num = 100, wave.height = 8, do.return =T) + 
  scale_fill_gradientn(colours = c("red","orange","yellow","green","blue","cyan","purple"))
ggsave(filename = "waddington.pdf", width = 10, height = 6)
waddingtonPlot <- function(layout, width = 1000, width.prop = 0.3, wave.num = 60, wave.height = 5, 
                           line.color = "black", line.type = 1, line.size = 1, line.alpha = 0.6, do.return = F){
  require(ggplot2)
  require(ggridges)
  # check layout is legal
  if(!(is.list(layout) &&
     sort(names(layout)) == sort(c("row", "row.scale", "col")) &&
     is.list(layout$col) &&
     is.vector(layout$row) &&
     is.vector(layout$row.scale) &&
     length(layout$row) == length(layout$col) && 
     length(layout$row) == 1+length(layout$row.scale) &&
     all(sapply(layout$col, length) == 2*layout$row))){
    stop("illegal layout!")
  }

  row.scale <- layout$row.scale
  waveform.nums <- layout$row
  waveform.curves <- lapply(1:length(waveform.nums), function(x){
    col.scale <- layout$col[[x]]
    sections <- unname(quantile(1:width, probs = col.scale/sum(col.scale), type = 1))
    curves <- sapply(1:(length(sections)/2), function(x) {
      cos.left <- cos((1:sections[2*x-1])*pi/sections[2*x-1])
      cos.right <- cos(pi+ (1:sections[2*x])*pi/sections[2*x])
      return(c(cos.left, cos.right))
    })
    return(unlist(curves)[1:width])
  })
  
  interpolate <- quantile(1:wave.num, probs = row.scale/sum(row.scale), type = 1)
  ggData <- NULL
  height_total <- sum(interpolate)
  ratio <- 1-width.prop
  for(i in 1:(length(waveform.nums)-1)){
    curve1 <- waveform.curves[[i]]
    curve2 <- waveform.curves[[i+1]]
    for(j in 0:(interpolate[i]-1)){
      alpha <- j/interpolate[i]
      curve <- (1-alpha)*curve1 + alpha*curve2 + 1
      height_j <- height_total - j - sign(i-1)*sum(interpolate[1:(i-1)])
      ggData <- rbind(ggData, data.frame(x = (1:width)*(1-ratio*height_j/height_total) + ratio*width*height_j/2/height_total, 
                                         y = unname(height_j), height = curve))
    }
  }
  p <- ggplot() + geom_ridgeline(data = ggData,
                            mapping = aes(x, y, height = height + 2, group = y, fill = y),
                            scale = wave.height, color = line.color, linetype = line.type, size = line.size, alpha = line.alpha) + theme_void()
  if(do.return){
    return(p)
  }else{
   print(p)
    return(NULL)
  }
}

