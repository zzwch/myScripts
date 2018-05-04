## 2d lightening
rotate <- function(x) t(apply(x, 2, rev))
n <- 10
library(grid)
mm <- tcrossprod(seq(1,0,length.out = n))
tmp1 <- sapply(col2rgb("orange")/255, function(x) 1-mm*(1-x))
tmp2 <- sapply(col2rgb("cyan")/255, function(x) 1-rotate(mm)*(1-x))
tmp3 <- sapply(col2rgb("purple")/255, function(x) 1-rotate(rotate(mm))*(1-x))
tmp4 <- sapply(col2rgb("grey")/255, function(x) 1-rotate(rotate(rotate(mm)))*(1-x))
colormat <- matrix(rgb(tmp1*tmp2*tmp3*tmp4),n)
grid.raster(colormat)

mapExpr2color <- function(gene1, gene2, colormat, scale = F){
  colormat <- apply(colormat, 2, rev)
  range1 <- quantile(gene1, probs = c(0.01, 0.99))
  if(!scale) range1[1] <- 0
  gene1[gene1 > range1[2]] <- range1[2]
  gene1[gene1 < range1[1]] <- range1[1]
  
  range2 <- quantile(gene2, probs = c(0.01, 0.99))
  if(!scale) range2[1] <- 0
  gene2[gene2 > range2[2]] <- range2[2]
  gene2[gene2 < range2[1]] <- range2[1]
  
  sapply(1:length(gene1), function(x) {
    colormat[round(1+(ncol(colormat)-1)*(gene2[x] - range2[1])/diff(range2)),
             round(1+(ncol(colormat)-1)*(gene1[x] - range1[1])/diff(range1))]
  })
}

gene1 <- "Kit"
gene2 <- "Runx1"
ggplot(data = data_annot) + 
  geom_point(mapping = aes(Main.wiPK44.tSNE_1, Main.wiPK44.tSNE_2), size = 2,
             color = mapExpr2color(as.numeric(data_expr[gene1, rownames(data_annot)]),
                                   as.numeric(data_expr[gene2, rownames(data_annot)]),
                                   colormat)) +
  ggtitle(label = paste0(gene1," vs ",gene2,"\n PCC = ", format(gcor[gene1,gene2], digits = 2), ", ",
                         "Pvalue = ", format(WGCNA::corPvalueStudent(gcor[gene1,gene2], nrow(annot)), 
                                             digits = 2, scientific = T)))
  
  grid.raster(colormat, x=4.5/5, y=1/4, width = 0.1, height = 0.1,interpolate = F)
  grid.text(x =4.5/5, y = 1/4-0.05, label = gene1, vjust = 1, gp = gpar(fontface = "italic"))
  grid.text(x =4.5/5-0.05, y = 1/4, label = gene2, vjust = -0.1, gp = gpar(fontface = "italic"), rot = 90)
