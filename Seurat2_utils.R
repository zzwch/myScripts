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
