writeSCENIC <- function(pbmc, file) {
  exprMat <- pbmc@data
  annot <- pbmc@meta.data[,"orig.ident",drop = F]
  save(exprMat, annot, file = file)
}

pairwise.topgene.corr <- function(dataMat, ngene = 500, method = "spearman"){
  order_gene <- apply(dataMat, 2, function(x){order(x, decreasing = T)[1:ngene]})
  resMat <- matrix(NA, nrow = ncol(dataMat), ncol = ncol(dataMat), dimnames = list(colnames(dataMat), colnames(dataMat)))
  for(i in 1:ncol(dataMat)){
    for(j in i:ncol(dataMat)){
      if(i == j){
        resMat[i,j] <- 0
      }else{
        data.x <- dataMat[,i,drop = T]
        data.y <- dataMat[,j,drop = T]
        sel_gene <- unique(c(order_gene[,i], order_gene[,j]))
        resMat[i,j] <- (1-cor(data.x[sel_gene], data.y[sel_gene], method = method))/2
        resMat[j,i] <- resMat[i,j]
      }
    }
  }
  return(resMat)
}

pairwise.topgene.equal <- function(dataMat, ngene = 500){
  order_gene <- apply(dataMat, 2, function(x){order(x, decreasing = T)[1:ngene]})
  resMat <- matrix(NA, nrow = ncol(dataMat), ncol = ncol(dataMat), dimnames = list(colnames(dataMat), colnames(dataMat)))
  for(i in 1:ncol(dataMat)){
    for(j in i:ncol(dataMat)){
      if(i == j){
        resMat[i,j] <- 0
      }else{
        equal <- length(intersect(order_gene[,i], order_gene[,j]))
        resMat[i,j] <- 1 - equal/(2*ngene-equal)
        resMat[j,i] <- resMat[i,j]
      }
    }
  }
  return(resMat)
}


focusTSNEPlot <- function(object, group.by = 'orig.ident', focus = NULL, each = F, colors.use = NULL,...){
  cats <- sort(unique(object@meta.data[[group.by]]))
  focus = intersect(focus, cats)
  if(each){
    sapply(focus, function(x){
      TSNEPlot(object, group.by =group.by, ...)
    })
  } else{
    if(is.null(names(colors.use)))
      grey_ind <- which(!(cats %in% focus))
    else
      grey_ind <- which(!(names(colors.use) %in% focus))
    colors.use[grey_ind] <- "grey"
    top <- object@meta.data[[group.by]] %in% focus
    toprow <- c(rownames(object@meta.data)[which(!top)],rownames(object@meta.data)[which(top)])
    object@dr$tsne@cell.embeddings <- object@dr$tsne@cell.embeddings[toprow,]
    TSNEPlot(object, group.by = group.by, colors.use = colors.use, ...)
  }
}

AddCCscore <- function(pbmc, th.g1s = 2, th.g2m = 2, org = c("mouse","human"), gene.max = NULL, add.dr = F, use.scale = F){
  geneset1 <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")
  geneset2 <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")
  org = org[1]
  if(org %in% c("human","hsa","hg19")){
    geneset1 <- toupper(geneset1)
    geneset2 <- toupper(geneset2)
  }
  use.data <- if(use.scale) pbmc@scale.data else pbmc@data
  use.data <- use.data[c(geneset1,geneset2),]
  use.data <- if(is.null(gene.max)) use.data else t(apply(use.data, 1, function(x) gene.max*x/max(x)))
  
  pbmc@meta.data$G1S.Score <- colMeans(use.data[geneset1, rownames(pbmc@meta.data)])
  pbmc@meta.data$G2M.Score <- colMeans(use.data[geneset2, rownames(pbmc@meta.data)])
  
  pbmc@meta.data$cc.phase[pbmc@meta.data$G1S.Score < th.g1s & pbmc@meta.data$G2M.Score < th.g2m] <- "Quiescent"
  pbmc@meta.data$cc.phase[!(pbmc@meta.data$G1S.Score < th.g1s & pbmc@meta.data$G2M.Score < th.g2m) & (pbmc@meta.data$G1S.Score < pbmc@meta.data$G2M.Score)] <- "G2/M"
  pbmc@meta.data$cc.phase[!(pbmc@meta.data$G1S.Score < th.g1s & pbmc@meta.data$G2M.Score < th.g2m) & (pbmc@meta.data$G1S.Score > pbmc@meta.data$G2M.Score) & pbmc@meta.data$G2M.Score < th.g2m] <- "G1"
  pbmc@meta.data$cc.phase[!(pbmc@meta.data$G1S.Score < th.g1s & pbmc@meta.data$G2M.Score < th.g2m) & (pbmc@meta.data$G1S.Score > pbmc@meta.data$G2M.Score) & pbmc@meta.data$G2M.Score > th.g2m] <- "S"
  pbmc@meta.data$cc.phase <- factor(pbmc@meta.data$cc.phase, levels = c("Quiescent","G1","S","G2/M"))
  
  
  if(add.dr){
    pbmc@dr$ccscore <- pbmc@dr$pca
    pbmc@dr$ccscore@cell.embeddings <- as.matrix(pbmc@meta.data[,c("G1S.Score","G2M.Score")])
  }
  pbmc
}
AddAVscore <- function(pbmc, org = c("mouse","human"), gene.max = NULL, add.dr = F, use.scale = F){
  geneset1 <- c("Efnb2","Dll4","Hey1","Gja4","Unc5b")
  geneset2 <- c("Ephb4","Nr2f2","Nrp2","Aplnr","Flt4")
  org = org[1]
  if(org %in% c("human","hsa","hg19")){
    geneset1 <- toupper(geneset1)
    geneset2 <- toupper(geneset2)
  }
  use.data <- if(use.scale) pbmc@scale.data else pbmc@data
  use.data <- use.data[c(geneset1,geneset2),]
  use.data <- if(is.null(gene.max)) use.data else t(apply(use.data, 1, function(x) gene.max*x/max(x)))
  
  pbmc@meta.data$Artery.Score <- colMeans(use.data[geneset1, rownames(pbmc@meta.data)])
  pbmc@meta.data$Venous.Score <- colMeans(use.data[geneset2, rownames(pbmc@meta.data)])
  if(add.dr){
    pbmc@dr$avscore <- pbmc@dr$pca
    pbmc@dr$avscore@cell.embeddings <- as.matrix(pbmc@meta.data[,c("Artery.Score","Venous.Score")])
  }
  pbmc
}
myGenePlot <- function(pbmc, gene1, gene2, use.scaled = F, group.by = NULL, cols.use = colors(),...){
  group.colors <- if(is.null(group.by)) pbmc@ident else pbmc@meta.data[[group.by]]
  plot(FetchData(pbmc, vars.all = gene1, use.scaled = use.scaled)[,1], FetchData(pbmc, vars.all = gene2,use.scaled = use.scaled)[,1],
       pch = 19, cex=1.6, col = cols.use[as.factor(group.colors)], xlab = gene1, ylab = gene2, ...)
}
myFeaturePlot <- function(pbmc, features.plot, nrow = NULL, ncol = NULL, dr = c("tsne","pca","ccscore","avscore"), cc.args = list(th.g1s = 2, th.g2m = 2),...){
  require(ggplot2)
  require(gridExtra)
  dr <- dr[1]
  ggData <- as.data.frame(cbind(pbmc@dr[[dr]]@cell.embeddings,FetchData(pbmc, features.plot)))
  colnames(ggData) <- c(colnames(pbmc@dr[[dr]]@cell.embeddings),gsub("-",".",features.plot))
  # print(feature.tmp)
  # ggData[,feature.tmp] <- t(pbmc@data[feature.tmp,])
  if(dr == 'tsne') {
    xx <- "tSNE_1"
    yy <- "tSNE_2"
  }
  if(dr == 'pca'){
    xx <- "PC1"
    yy <- "PC2"
  }
  if(dr == 'ccscore'){
    xx <- "G1S.Score"
    yy <- "G2M.Score"
  }
  if(dr == 'avscore'){
    xx <- "Artery.Score"
    yy <- "Venous.Score"
  }
  ggl <- lapply(features.plot, function(feature){
    p <- ggplot(ggData) + geom_point(mapping = aes_string(x = xx, y = yy, color = gsub("-",".",feature)), size = 2) + 
      scale_color_gradientn(colours = c("grey","yellow","red")) + 
      xlab(label = xx) + ylab(label = yy) +
      theme(legend.title = element_blank()) + ggtitle(feature) 
    if(dr == "ccscore"){
      th.g1s <- cc.args$th.g1s
      th.g2m <- cc.args$th.g2m
      ccx <- ceiling(max(pbmc@meta.data$G1S.Score))
      ccy <- ceiling(max(pbmc@meta.data$G2M.Score))
      p <- p + geom_linerange(mapping = aes_(x = th.g1s, ymin = 0, ymax = th.g1s)) + 
        geom_segment(mapping = aes_(x = 0, y = th.g2m, xend = ccx, yend = th.g2m)) + 
        geom_segment(mapping = aes_(x = th.g1s, y = th.g2m, xend = ccx, yend = ccy)) +
        geom_text(mapping = aes_(th.g1s, quote(0.3), label = quote("Quiescent"), hjust = 1.1)) +
        geom_text(mapping = aes_(ccx, quote(0.3), label = quote("G1"), hjust = 1.1)) +
        geom_text(mapping = aes_(ccx, ccy-2, label = quote("S"), hjust = 1.1)) +
        geom_text(mapping = aes_(th.g1s/2, ccy-1, label = quote("G2M"), hjust = 1.1))
    }
    p
  })
  grid.arrange(grobs = ggl, nrow= nrow,ncol = ncol)
}

myCCPlot <- function(pbmc, group.by = NULL, th.g1s = 2, th.g2m = 2, text.size = 3, cols.use = palette(),...){
  group.colors <- sort(unique(if(is.null(group.by)) pbmc@ident else pbmc@meta.data[[group.by]]))
  require(dplyr)
  require(ggplot2)
  ccData <- pbmc@meta.data 
  ccx <- ceiling(max(ccData$G1S.Score))
  ccy <- ceiling(max(ccData$G2M.Score))
  p1 <- ggplot(ccData) + geom_point(mapping = aes_string("G1S.Score", "G2M.Score", color = group.by), size = 3) + 
    xlab(label = "G1/S phase score") + ylab("G2/M phase score") + ggtitle("Cell Cycle Analysis") +
    scale_color_manual(values = cols.use[as.factor(group.colors)]) +
    #scale_x_continuous(expand = c(0, 0), limits = c(0, ccx)) + 
    #scale_y_continuous(expand = c(0, 0), limits = c(0, ccy)) +
    geom_linerange(mapping = aes_(x = th.g1s, ymin = 0, ymax = th.g1s)) + 
    geom_segment(mapping = aes_(x = 0, y = th.g2m, xend = ccx, yend = th.g2m)) + 
    geom_segment(mapping = aes_(x = th.g1s, y = th.g2m, xend = ccx, yend = ccy)) +
    geom_text(mapping = aes_(th.g1s, quote(0.3), label = quote("Quiescent"), hjust = 1.1)) +
    geom_text(mapping = aes_(ccx, quote(0.3), label = quote("G1"), hjust = 1.1)) +
    geom_text(mapping = aes_(ccx, ccy-2, label = quote("S"), hjust = 1.1)) +
    geom_text(mapping = aes_(th.g1s/2, ccy-1, label = quote("G2M"), hjust = 1.1)) +
    theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  
  ccData <- pbmc@meta.data %>% count_(vars = c(group.by, "cc.phase")) %>% group_by_(group.by) %>% arrange(desc(cc.phase)) %>%
    mutate(pct = n/sum(n),
           ypos = cumsum(pct) - 0.5*pct)
  p2 <- ggplot(ccData) + 
    geom_bar(mapping = aes_string(x = group.by, y = "pct", fill = "cc.phase"), stat = "identity", width = 0.8) + 
    geom_text(mapping = aes_string(x = group.by, y= "ypos",label = quote(paste0(ccData$n,", ",sprintf("%1.1f", 100*ccData$pct),"%"))), size = text.size) +
    ggtitle(label = "Cell cycle Distribution of Cells") + theme(plot.title = element_text(hjust = 0.5))
#  pheatmap::pheatmap(expr, cluster_cols = F, cluster_rows = F, annotation_col = annot, border_color = NA)
  require(gridExtra)
  grid.arrange(grobs = list(p1,p2), ncol = 2)
}
                                                          
bpGOEnrich <- function(pbmc.markers){
  require(BiocParallel)
  
  bplapply(sort(unique(pbmc.markers$cluster)), FUN = function(i){
    require(clusterProfiler)
    require(org.Hs.eg.db)
    require(stringr)
    require(ggplot2)
    gogenes <- subset(pbmc.markers, cluster == i)$gene
    
    engo <- enrichGO(gene         = gogenes,
                     OrgDb         = org.Hs.eg.db,
                     keytype       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
    png(filename = paste0("GO_BP results of Cluster",i,".png"), width = 600, height = 350)
    print(dotplot(engo, title = paste0("GO_BP results of Cluster",i)) + scale_y_discrete(labels = function(x) str_wrap(x, width = 60)))
    dev.off()
  }, BPPARAM = param)
  
}
myGOEnrich <- function(pbmc.markers, clusters = NULL, width = 600, height = 350){
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(stringr)
  if(is.null(clusters)){
    clusters <- sort(unique(pbmc.markers$cluster))
  }
  for(i in clusters){
    gogenes <- subset(pbmc.markers, cluster == i)$gene
    
    engo <- enrichGO(gene         = gogenes,
                     OrgDb         = org.Hs.eg.db,
                     keytype       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
    png(filename = paste0("GO_BP results of Cluster",i,".png"), width = width, height = height)
    print(dotplot(engo, title = paste0("GO_BP results of Cluster",i))+ scale_y_discrete(labels = function(x) str_wrap(x, width = 60)))
    dev.off()
  }
}
# convert seurat2 Object to Rdata needed by sciDV (a in-house shiny-server) 
seurat2idv <- function(pbmc = NULL, subcat = NULL, colors = NULL, save.file = NULL){
  data_expr <- pbmc@data
  data_annot <- pbmc@meta.data
  data_subcat <- intersect(subcat, colnames(pbmc@meta.data))
  data_color <- colors
  data_coord <- list()
  if(!is.null(pbmc@dr$pca)){
    data_annot <- cbind(data_annot, pbmc@dr$pca@cell.embeddings[,c(1,2)])
    data_coord <- c(data_coord, list(PCA = colnames(pbmc@dr$pca@cell.embeddings[,c(1,2)])))
  }
  if(!is.null(pbmc@dr$tsne)){
    data_annot <- cbind(data_annot, pbmc@dr$tsne@cell.embeddings[,c(1,2)])
    data_coord <- c(data_coord, list(tSNE = colnames(pbmc@dr$tsne@cell.embeddings[,c(1,2)])))
  }
  
  save(data_expr, data_annot, data_subcat, data_color, data_coord, file = save.file)
  return()
}
