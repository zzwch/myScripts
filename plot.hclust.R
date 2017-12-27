options(stringsAsFactors = F)
hc <- hclust(dist(USArrests), "ave")
plot(hc)
plot(hc, hang = -1)

if(!require(ggdendro)){
  install.packages("ggdendro")
  library(ggdendro)
}

hc_data <- dendro_data(hc, type = "triangle")
str(hc_data)
# List of 4
# $ segments   :'data.frame':	98 obs. of  4 variables:
#   ..$ x   : num [1:98] 17.99 4.77 1.5 1.5 4.77 ...
# ..$ y   : num [1:98] 152.3 77.6 38.5 38.5 77.6 ...
# ..$ xend: num [1:98] 4.77 1.5 1 2 8.03 ...
# ..$ yend: num [1:98] 77.6 38.5 0 0 44.3 ...
# $ labels     :'data.frame':	50 obs. of  3 variables:
#   ..$ x    : num [1:50] 1 2 3 4 5 6 7 8 9 10 ...
# ..$ y    : num [1:50] 0 0 0 0 0 0 0 0 0 0 ...
# ..$ label: Factor w/ 50 levels "Florida","North Carolina",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ leaf_labels: NULL
# $ class      : chr "hclust"
# - attr(*, "class")= chr "dendro"

hc_data$labels$color <- sample(colors(), length(hc_data$labels$label))
library(ggplot2)
ggplot() + 
  geom_segment(data = hc_data$segments, 
               mapping = aes(x, y, xend = xend, yend = yend)) +
  geom_text(data = hc_data$labels, 
             mapping = aes(x, y-1, label = label, color = color), 
             angle = 45, lineheight = 0, hjust = 1) +
  xlim(-3,NA) + ylim(-40,NA) +
  theme_dendro() + 
  theme(legend.position = "none")


library(ggtree)
library(ape)
hc_phylo <- as.phylo(hc)
str(hc_phylo)
ggtree(hc_phylo,color = "purple", layout = "fan", open.angle = 100,
       size = 0.5, linetype = "dotted", ladderize = T,
       branch.length = "branch.length") +
  geom_tippoint(color = "magenta", shape = "@", size =1) +
  geom_tiplab2(aes(angle = angle), size = 2, color = "orange") +
  theme_tree(bgcolor = "black")
