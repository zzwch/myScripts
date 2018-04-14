options(stringsAsFactors = F)
rm(list = ls())

getAmplifiedLibrary <- function(content = 1e4, cycles = 12, amplification_prob = 1){
  library <- 1:content
  # for(cy in 1:cycles){
  #   content_ <- sapply(X = amplification_prob,FUN = function(x){
  #     sample(c(T,F), 1, prob = c(x,1-x))
  #   })
  # }
  # 
  # distribution_fun <- get(paste0("r",distribution))
  # dup_rate_per_region <- distribution_fun(length(non_redundant_regions))
  # dup_rate_per_region <- round(100*(dup_rate_per_region - min(dup_rate_per_region)))
  if(all(amplification_prob >= 1))
    return(rep(library, each = 1*2^cycles))
}


# variable cell_content and variable sequencing_depth
cell_content <- 10^seq(4, 7, 0.2)
xx <- list()
yy <- list()
for(i in 1:length(cell_content)){
  library<- getAmplifiedLibrary(content = cell_content[i], amplification = 1, cycle = 5)
  
  xx[[i]] <- 10^seq(1,5,0.2)
  yy[[i]] <- xx[[i]]*0
  for(j in 1:length(xx[[i]])){
    depth <- xx[[i]][j]
    seq <- sample(library, depth)
    yy[[i]][j] <- sum(duplicated(seq))/depth
  }
}
ggData <- data.frame(seq_depth = unlist(xx), dup_rate = unlist(yy), cell_content = rep(cell_content, each = length(xx[[1]])))

library(ggplot2)
ggplot(ggData) + geom_line(mapping = aes(seq_depth, 100*dup_rate, group = cell_content, color = log10(cell_content)), size = 1) +
  scale_color_gradientn(colours = rev(c("red","orange","yellow","green","cyan","blue","purple"))) +
  labs(color = "DNA content,\nlog10 scaled")+
  xlab("Sequencing depth (Reads)") + ylab("Duplicated rate (%)") + 
  ggtitle(paste0("Changes of Duplication rates \nunder different Sequencing depth and DNA content")) + 
  geom_text(mapping = aes(0,100, label = "Assumption: \nDuring PCR cycles, reads are amplificated \nlinearly & completely."), 
            size = 4, color = "red", hjust = 0, vjust = 1)+
  ylim(c(0,100)) + theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_text(face = "bold"), 
        title = element_text(face="bold"), text = element_text(colour = "darkblue"))
