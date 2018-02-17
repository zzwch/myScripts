waddingtonPlot(valleys = c(1,1,2,4,8), 
               valleys_layout = lapply(c(1,1,2,4,8), function(x) {rep(1,2*x)}),
               interpolator = c(2,5,8,20) *2,
               line.type = 2, 
               line.size = 0.1, 
               top_scale = 0.25,
               theme.void = T, do.return = T) + 
  annotation_custom(grob = grid::rasterGrob(image = EBImage::readImage(files = "微信图片_20180217191752.jpg")),
                    450, 550, 77, 85)
ggsave(filename = "waddington.toy.pdf", width = 8, height = 6)
########### Main Function Body #############
waddingtonPlot <- function(valleys = c(1,2,4), 
                           valleys_layout = list(c(1,1),               # valley = 1
                                                 c(1,1, 1,1),         # valley = 2
                                                 c(1,1, 1,1, 1,1, 1,1) # valley = 4
                                                ),
                           interpolator = c(20,20), 
                           valley_width = 1000,
                           top_scale = 0.3,
                           ridge.height = 5, 
                           ridge.colors = c("red","orange","yellow","green","blue","cyan","purple"),
                           line.color = "black", line.type = 1, line.size = 1, line.alpha = 0.2, 
                           theme.void = T, hide.legend = T, do.return = F){
  require(ggplot2)
  require(ggridges)

  waveform.curves <- lapply(1:length(valleys), function(x){
    bias <- valleys_layout[[x]]
    sections <- unname(quantile(1:valley_width, probs = bias/sum(bias), type = 1))
    curves <- sapply(1:(length(sections)/2), function(x) {
      cos.left <- cos((1:sections[2*x-1])*pi/sections[2*x-1])
      cos.right <- cos(pi+ (1:sections[2*x])*pi/sections[2*x])
      return(c(cos.left, cos.right))
    })
    return(unlist(curves)[1:valley_width])
  }) # get the main wave curves of each one in valleys
  
  ggData <- NULL
  waves_sum <- sum(interpolator)
  ratio <- 1-top_scale
  for(i in 1:(length(valleys)-1)){
    curve1 <- waveform.curves[[i]]
    curve2 <- waveform.curves[[i+1]]
    for(j in 0:(interpolator[i]-1)){
      alpha <- j/interpolator[i]
      curve <- (1-alpha)*curve1 + alpha*curve2 + 1
      waves_cur <- waves_sum - j - sign(i-1)*sum(interpolator[1:(i-1)])
      ggData <- rbind(ggData, data.frame(x = (1:valley_width)*(1-ratio*waves_cur/waves_sum) + ratio*valley_width*waves_cur/2/waves_sum, 
                                         y = unname(waves_cur), shift = curve))
    }
  }
  p <- ggplot() + 
    geom_ridgeline(data = ggData,
                   mapping = aes(x, y, height = shift + 2, group = y, fill = y),
                   scale = ridge.height,
                   color = line.color,
                   linetype = line.type, 
                   size = line.size, 
                   alpha = line.alpha) +
    scale_fill_gradientn(colours = ridge.colors)
  if(theme.void)  p <- p + theme_void()
  if(hide.legend) p <- p + theme(legend.position = "none")
  
  
  if(do.return){
    return(p)
  }else{
    print(p)
    return(NULL)
  }
}

