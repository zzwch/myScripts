pal_dolphin <- c('pink','#FF00AE','#A020F1','#000000','#0403E5','#FF8C01','#8B0101','#007502','#FE0000','#FFFF01','#FF99CB','#4A95FB','#61FE69','#9A7A01','#017F8B','#05FDFF','grey','grey','wheat','orange4','red','magenta','gold','#66A61E','skyblue','#D0B38A','#A3D171')
pal_colorbrewer_spectral <- RColorBrewer::brewer.pal(11, "Spectral")
pal_colorbrewer_paired <- RColorBrewer::brewer.pal(12, "Paired")

specScore <- function(){
  # to be done!!!!!!
  # internal function to calculate AUC values
  AUCMarkerTest <- function(data1, data2, mygenes, print.bar = TRUE) {
    myAUC <- unlist(x = lapply(
      X = mygenes,
      FUN = function(x) {
        return(DifferentialAUC(
          x = as.numeric(x = data1[x, ]),
          y = as.numeric(x = data2[x, ])
        ))
      }
    ))
    myAUC[is.na(x = myAUC)] <- 0
    if (print.bar) {
      iterate.fxn <- pblapply
    } else {
      iterate.fxn <- lapply
    }
    avg_diff <- unlist(x = iterate.fxn(
      X = mygenes,
      FUN = function(x) {
        return(
          ExpMean(
            x = as.numeric(x = data1[x, ])
          ) - ExpMean(
            x = as.numeric(x = data2[x, ])
          )
        )
      }
    ))
    toRet <- data.frame(cbind(myAUC, avg_diff), row.names = mygenes)
    toRet <- toRet[rev(x = order(toRet$myAUC)), ]
    return(toRet)
  }
}
# internal function to calculate spec values
DifferentialPrediction <- function(x, y, measure = "spec") {
  prediction.use <- ROCR::prediction(
    predictions = c(x, y),
    labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))),
    label.ordering = 0:1
  )
  perf.use <- ROCR::performance(prediction.obj = prediction.use, measure = measure)
  return(round(x = max(perf.use@y.values[[1]]), digits = 3))
}

autofilterBarcode <- function(umi, plot = T, do.return = F, p.cutoff = 0.5){
  # do.return - return ggplot object or 'not' (return barcodes filtered umi matrix)
  # plot - plot the ggplot object or not
  require(mixtools)
  require(ggplot2)
  mixmdl <- normalmixEM(log10(colSums(umi)+1), k = 2)
  ind_valid <- which(mixmdl$posterior[,which.max(mixmdl$mu)] > 1-p.cutoff)
  new_umi <- umi[, ind_valid]
  
  plot_mix_comps <- function(x, mu, sigma, lam) {
    lam * dnorm(x, mu, sigma) * length(x)
  }
  curve_colors <- c("grey","red")
  if(mixmdl$mu[1] > mixmdl$mu[2]) curve_colors <- c("red","grey")
  p1 <- ggplot(data.frame(x = mixmdl$x)) +
    geom_histogram(aes(x, ..count..), bins = 100, colour = "black", 
                   fill = "white") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                  colour = curve_colors[1], lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                  colour = curve_colors[2], lwd = 1.5) +
    labs(x = "nUMI (log10 scaled)", y = "nCells", 
         title = paste0("Histogram of nUMI \n",
                        "Valid nUMI = ",sum(new_umi),", ", round(sum(new_umi)/sum(umi), digits = 4) * 100,"%;\n",
                        "Valid nCell = ",ncol(new_umi),"; Average nUMI = ",round(sum(new_umi)/ncol(new_umi), digits = 1))) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  ggData2 <- cbind(reshape2::colsplit(colnames(umi), "_sc", names = c("Library", "Barcode")), valid = "0")
  ggData2$valid[ind_valid] <- "1"
  p2 <- ggplot(ggData2)+ geom_tile(mapping = aes(x = Barcode, y = Library, fill = valid)) +
    scale_fill_manual(values = c("grey","red"))+labs(title = "Barcode Usage") + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  if(plot) gridExtra::grid.arrange(p1,p2, nrow = 1)
  if(do.return) return(list(p1,p2))
  
  return(new_umi)
}

dptPlot <- function(dpt, root = NULL, paths_to = integer(0L), dcs = 1:2, 
                    divide = integer(0L), w_width = 0.1, annot, col_by = "dpt", shape_by = "branch",
                    col_path = c("red","darkblue"), col_tip = "red", 
                    col = NULL, legend_main = col_by){
  require(ggplot2)
  dpt_flat <- branch_divide(dpt, divide)
  if (!is.null(root) && length(root) < 1L) 
    stop("root needs to be specified")
  root <- {if (is.null(root)) min(dpt_flat@branch[, 1], na.rm = TRUE)
           else as.integer(root)}
  paths_to <- as.integer(paths_to)
  if (length(root) > 1L && length(paths_to) > 0L) 
    stop("(length(root), length(paths_to)) needs to be (1, 0-n) or (2-n, 0), but is (", 
         length(root), ", ", length(paths_to), ")")
  stopifnot(length(dcs) %in% 2:3)
  if (length(root) > 1L && length(paths_to) == 0L) {
    paths_to <- root[-1]
    root <- root[[1]]
  }
  
  # xy
  eiv <- dpt_flat@dm@eigenvectors[,dcs,drop = F]
  #time
  branch_idx <- dpt_flat@branch[, 1L] == root
  stopifnot(any(branch_idx, na.rm = T))
  tip_cells <- which(branch_idx & dpt_flat@tips[, 1L])
  if (length(tip_cells) == 0L) tip_cells <- which(branch_idx)
  pt <- dpt[tip_cells[[1L]], ]
  #branch
  brc <- dpt_flat@branch[,1]
  
  # annot
  annot <- pData(dpt_flat@dm@data_env$data)
  ggData <- as.data.frame(cbind(eiv, branch = paste0("branch",brc), dpt = pt, annot))
  ggplot(ggData) + geom_point(mapping = aes_string(x = "DC1", y = "DC2", color = col_by, shape = shape_by), size = 2) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}
seuratPCAtopGene <- function(pbmc, pcs = 1, top = 20, pc.use.full = T, plot = T, plot.scale = F, row.scale = F){
  pca <- pbmc@dr$pca
  pc <- {if(pc.use.full) pca@gene.loadings.full else pca@gene.loadings}
  phData <- {if(plot.scale) pbmc@scale.data else pbmc@data}
  genescores <- sort(pc[,pcs], decreasing = (top > 0))[1:abs(top)]
  
  pheatmap::pheatmap(phData[names(genescores), order(pca@cell.embeddings[,pcs]), drop = F], cluster_rows = F, cluster_cols = F,
                     scale = ifelse(row.scale, "row","none"),
                     annotation_col = as.data.frame(pca@cell.embeddings[,pcs, drop = F]), show_colnames = F)
  return(genescores)
}
ggGeneData <- function(expr, annot, genes, annot.col, annot.col.id.vars = NULL){
  require(reshape2)
  annot.col = unique(c(annot.col, annot.col.id.vars))
  if(length(intersect(annot.col, colnames(annot))) == 0) stop("none col selected")
  samples <- intersect(rownames(annot), colnames(expr))
  annot <- annot[samples, annot.col, drop = F]
  expr <- expr[genes,samples, drop = F]
  if(is.null(annot.col.id.vars)){
    ggData <- melt(as.data.frame(cbind(annot, t(expr))))
  }else{
    ggData <- melt(as.data.frame(cbind(annot, t(expr))), id.vars = annot.col.id.vars)
  }
}
writeList <- function(l, file = "clipboard", sep = "\t"){
  con <- file(file, open = "w")
  lapply(1:length(l), function(i){
    if(!is.null(names(l)))
      writeLines(names(l)[i], con = con, sep = sep)
    writeLines(l[[i]], con = con, sep = sep)
    writeLines("", con = con, sep = "\n")
  })
  close(con)
}
myRowScale <- function(expr, cap = 10, quantile.probs = 0.99){
  t(apply(expr, 1, function(x){
    robustMax <- quantile(x, probs = quantile.probs)
    x[x > robustMax] <- robustMax
    x <- cap*x/robustMax
  }))
}
capping <- function(m, max = 2, min =-2){
  m[m > max] <- max
  m[m < min] <- min
  m
}
matScale <- function(mat, scale.by = c("row","col")){
  res <- mat
  scale.by <- scale.by[1]
  res[,] <- if(scale.by == "row") t(scale(t(mat))) else scale(mat)
  return(res)
}
rowSmooth <- function(mat, n = NULL, proportion = 0.2){
  if(is.na(n) || is.null(n)){
    n <- ceiling(ncol(mat) * proportion)
  }
  if(n == 1) stop("nothing should be smoothed!")
  res <- mat
  for(i in 1:nrow(mat)){
    res[i,] <- vectorSmooth(as.numeric(mat[i,,drop = T]), n)
  }
  return(res)
}
vectorSmooth <- function(x, n){
  sapply(1:length(x), function(i){
    m <- 1:n + i - ceiling(n/2)
    m <- m[m > 0 & m <= length(x)]
    mean(x[m])
  })
}

reduceColor <- function(colors, annot){
  col = intersect(names(colors), colnames(annot))
  for(i in col){
    if(is.factor(annot[[i]])){
      colors[[i]] <- colors[[i]][levels(annot[[i]])]
      
    }else{
      colors[[i]] <- colors[[i]][unique(sort(annot[[i]]))]
    }
  }
  return(colors)
}
writeSCENIC <- function(pbmc, annot_col,file) {
  exprMat <- pbmc@data
  annot <- pbmc@meta.data[,annot_col,drop = F]
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
gradient2d <- function(fourcolor = c(lefttop = "orange", righttop = "cyan", leftbottom = "grey", rightbottom = "purple"), n = 10, plot = F){
  require(grid)
  if(is.null(names(fourcolor))) names(fourcolor) <- c("lefttop", "righttop","leftbottom", "rightbottom") 
  rotate <- function(x) t(apply(x, 2, rev))
  mm <- tcrossprod(seq(1,0,length.out = n))
  tmp1 <- sapply(col2rgb(fourcolor[["lefttop"]])/255, function(x) 1-mm*(1-x))
  tmp2 <- sapply(col2rgb(fourcolor[["righttop"]])/255, function(x) 1-rotate(mm)*(1-x))
  tmp3 <- sapply(col2rgb(fourcolor[["rightbottom"]])/255, function(x) 1-rotate(rotate(mm))*(1-x))
  tmp4 <- sapply(col2rgb(fourcolor[["leftbottom"]])/255, function(x) 1-rotate(rotate(rotate(mm)))*(1-x))
  
  if(!plot) colormat <- matrix(rgb(tmp1*tmp2*tmp3*tmp4),n)
  else{
    plot.new()
    grid.raster(colormat)
  }
}
lightening2d <- function(data_annot, axis.x, axis.y, data_expr, gene1, gene2, pt.size = 2, cor.title = F,
                         fourcolor = c(lefttop = "orange", righttop = "cyan", leftbottom = "grey", rightbottom = "purple"), n = 10,
                         x.2d = 0.9, y.2d = 0.25, w.2d = 0.1, h.2d = 0.1, interpolate.2d = F){
  require(ggplot2)
  require(WGCNA)
  require(Seurat)
  ## 2d lightening
  colormat <- gradient2d(fourcolor, n = n, plot = F)
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

  p <- ggplot(data = data_annot) + #theme_classic()+
    geom_point(mapping = aes_string(axis.x, axis.y), size = pt.size,
               color = mapExpr2color(as.numeric(data_expr[gene1, rownames(data_annot)]),
                                     as.numeric(data_expr[gene2, rownames(data_annot)]),
                                     colormat))
  if(cor.title){
    gcor <- WGCNA::cor(t(data_expr[c(gene1,gene2), ]))
    p <- p + ggtitle(label = paste0(gene1," vs ",gene2,"\n PCC = ", format(gcor[gene1,gene2], digits = 2), ", ",
                                    "Pvalue = ", format(WGCNA::corPvalueStudent(gcor[gene1,gene2], nrow(data_annot)), 
                                                        digits = 2, scientific = T)))
  }else{
    p <- p + ggtitle(label = paste0(gene1," vs ",gene2))
  }
  print(p)
  grid.raster(colormat, x= x.2d, y=y.2d, width = w.2d, height = h.2d,interpolate = interpolate.2d)
  grid.text(x =x.2d, y =y.2d-h.2d/2, label = gene1, vjust = 1, gp = gpar(fontface = "italic"))
  grid.text(x =x.2d-w.2d/2, y = y.2d, label = gene2, vjust = -0.2, gp = gpar(fontface = "italic"), rot = 90)
}

degcomp <- function(marker1, marker2){
  marker1 <- subset(marker1, p_val < 0.05)
  marker1.up <- subset(marker1, avg_logFC > 0)
  marker1.dn <- subset(marker1, avg_logFC < 0)
  
  marker2 <- subset(marker2, p_val < 0.05)
  marker2.up <- subset(marker2, avg_logFC > 0)
  marker2.dn <- subset(marker2, avg_logFC < 0)
  
  m <- matrix(NA, 2, 3)
  rownames(m) <- c("up", "dn")
  colnames(m) <- c("marker1", "marker2", "intersect")
  m[1,1] <- nrow(marker1.up)
  m[1,2] <- nrow(marker2.up)
  m[1,3] <- length(intersect(rownames(marker1.up), rownames(marker2.up)))
  m[2,1] <- nrow(marker1.dn)
  m[2,2] <- nrow(marker2.dn)
  m[2,3] <- length(intersect(rownames(marker1.dn), rownames(marker2.dn)))
  return(m)
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
AddAVscore <- function(pbmc, org = c("mouse","human"), gene.max = NULL, add.dr = F, use.scale = F, 
                       gs1 = c("Efnb2","Dll4","Hey1","Gja4","Unc5b"), gs2 = c("Ephb4","Nr2f2","Nrp2","Aplnr","Flt4")){
  geneset1 <- gs1
  geneset2 <- gs2
  
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
seuratAddDPT <- function(pbmc, dpt){
  pbmc@dr$dpt <- pbmc@dr$tsne
  tmp <- as.matrix(dpt@dm@eigenvectors)
  rownames(tmp) <- rownames(pData(dpt@dm@data_env$data))
  pbmc@dr$dpt@cell.embeddings <- tmp[rownames(pbmc@dr$tsne@cell.embeddings),]
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
  if(dr == 'dpt'){
    xx <- "DC1"
    yy <- "DC2"
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

ggHsegments <- function(tab, y = 0, group){
  tab <- tab[tab > 0] # tab must be 1-d and derived from table()
  seg <- matrix(NA, length(tab), 5, dimnames = list(as.character(names(tab)), c(group,"x","y","xend","yend")))
  seg <- as.data.frame(seg)
  seg[,group] <- rownames(seg)
  seg[,"x"] <- c(0,cumsum(tab)[-length(tab)]) + 0.5
  seg[,"xend"] <- cumsum(tab) +0.5
  seg[,"y"] <- y
  seg[,"yend"] <- y
  return(seg)
}
myBarPlot <- function(pbmc, features.plot, title = NULL, nrow = NULL, ncol = NULL, group.by = NULL, bar.width = 0.5, x.order = NULL,
                      ymax = c("auto",10), facet.scales = c("fixed", "free","free_x","free_y"),cols.use = pal_dolphin, use.scale = F){
  require(ggplot2)
  require(reshape2)
  if(is.null(group.by)) {
    pbmc@meta.data$ident <- pbmc@ident
    group.by <- "ident"
  }
  facet.scales <- facet.scales[1]
  expr <- if(use.scale) pbmc@scale.data else pbmc@data
  ggData <- cbind(pbmc@meta.data[,group.by, drop = F], t(expr)[rownames(pbmc@meta.data), features.plot, drop = F])
  ggData$my.rownames <- rownames(ggData)
  set.seed(666)
  ggData$my.order <- {if(is.null(x.order)) rank(ggData[[group.by]],ties.method = "random") else x.order}
  ggData <- melt(ggData, id.vars = c("my.rownames","my.order", group.by), variable.name = "Select_feature")
  ggSeg <- ggHsegments(table(pbmc@meta.data[[group.by]]), y = 0, group = group.by)
  
  p <- ggplot(ggData) + 
    geom_bar(mapping = aes_string(x = "my.order",y = "value", fill = group.by), width = bar.width, stat = "identity") +
    geom_segment(data = ggSeg, mapping = aes_string(x = "x", y = "y",xend = "xend", yend = "yend",color = group.by)) + 
    facet_wrap(~Select_feature, nrow = nrow, ncol = ncol, scales = facet.scales) +
    labs(x = NULL, y = NULL, title = title) + 
    scale_fill_manual(values = cols.use) + scale_color_manual(values = cols.use) +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = NA, colour = NA),
          strip.placement = "inside",
          strip.text = element_text(face = "bold.italic"))
  
  ymax <- ymax[1]
  if(ymax == "auto"){
    #p + coord_cartesian(xlim = c(-1,nrow(pbmc@meta.data)+1), ylim = c(0,ceiling(max(ggData$value))),expand = F)
    p + coord_cartesian(xlim = c(-1,nrow(pbmc@meta.data)+1),expand = F)
  }else{
    p + coord_cartesian(xlim = c(-1,nrow(pbmc@meta.data)+1), ylim = c(0,ymax),expand = F)
  }
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
    geom_linerange(mapping = aes_(x = th.g1s, ymin = 0, ymax = th.g2m)) + 
    geom_segment(mapping = aes_(x = 0, y = th.g2m, xend = ccx, yend = th.g2m)) + 
    geom_segment(mapping = aes_(x = th.g1s, y = th.g2m, xend = ccx, yend = ccy)) +
    geom_text(mapping = aes_(th.g1s, quote(0.3), label = quote("Quiescent"), hjust = 1.1)) +
    geom_text(mapping = aes_(ccx, quote(0.3), label = quote("G1"), hjust = 1.1)) +
    geom_text(mapping = aes_(ccx, ccy-2, label = quote("S"), hjust = 1.1)) +
    geom_text(mapping = aes_(th.g1s/2, ccy-1, label = quote("G2M"), hjust = 1.1)) +
    theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  
  ccData <- pbmc@meta.data %>% count_(vars = c(group.by, "cc.phase")) %>% group_by_(group.by) %>% dplyr::arrange(desc(cc.phase)) %>%
    dplyr::mutate(pct = n/sum(n),
           ypos = cumsum(pct) - 0.5*pct)
  p2 <- ggplot(ccData) + 
    geom_bar(mapping = aes_string(x = group.by, y = "pct", fill = "cc.phase"), stat = "identity", width = 0.8) + 
    geom_text(mapping = aes_string(x = group.by, y= "ypos",label = quote(paste0(ccData$n,", ",sprintf("%1.1f", 100*ccData$pct),"%"))), size = text.size) +
    ggtitle(label = "Cell cycle Distribution of Cells") + 
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust =1, vjust = 1))
  #  pheatmap::pheatmap(expr, cluster_cols = F, cluster_rows = F, annotation_col = annot, border_color = NA)
  require(gridExtra)
  grid.arrange(grobs = list(p1,p2), ncol = 2)
}
myAVPlot <- function(pbmc, group.by = NULL, continuous = F, gene.max = 10, use.lim = T, pt.size = 2, text.size = 3, cols.use = palette(), aggregate.fun = mean, ...){
  # atery venous score
  require(ggplot2)
  require(gridExtra)
  library(ggrepel)
  pbmc@meta.data$ident <- pbmc@ident[rownames(pbmc@meta.data)]
  if(is.null(group.by)) group.by <- "ident"
  p1 <- ggplot(data = pbmc@meta.data, mapping = aes_string(x = "Venous.Score", y = "Artery.Score", color = group.by)) + 
    geom_point(size = pt.size) + 
    coord_fixed(ratio = 1)
  if(use.lim) p1 <- p1 + xlim(0,gene.max) + ylim(0,gene.max)
  if(continuous) p1 + scale_color_gradientn(colours = cols.use)
  else{
    p1 <- p1 + scale_color_manual(values = cols.use)
    
    aveArtery <- aggregate(pbmc@meta.data$Artery.Score, list(Average = pbmc@meta.data[[group.by]]), function(y){aggregate.fun(y)})[,2]
    aveVenous <- aggregate(pbmc@meta.data$Venous.Score, list(Average = pbmc@meta.data[[group.by]]), function(y){aggregate.fun(y)})[,2]
    
    av_ggData <- data.frame(Venous_average = aveVenous, Artery_average = aveArtery, Average = as.character(unique(sort(pbmc@meta.data[[group.by]]))))
    p2 <- ggplot(av_ggData,aes(x = Venous_average, y = Artery_average)) + 
      geom_point(aes(color = Average), size = pt.size) + 
      geom_label_repel(aes(label = Average, color = Average)) + 
      scale_color_manual(values = cols.use) + 
      coord_fixed(ratio = 1) + guides(color = guide_legend(title = group.by))
    if(use.lim) p2 <- p2 + xlim(0,gene.max) + ylim(0,gene.max)
    grid.arrange(p1,p2, nrow = 1)
  } 
}
myHeatmap <- function(object, use.scaled = T, genes.use = NULL, annot.use = NULL, annot.colors = NULL,
                      disp.min = -2.5, disp.max = 2.5, pdf.file = NULL, pdf.width, pdf.height ){
  if(use.scaled){
    data.use <- object@scale.data
  }else{
    data.use <- object@data
  }
  
  data.use <- data.use[intersect(genes.use, rownames(data.use)), order(object@ident)]
  data.use[data.use > disp.max] <- disp.max
  data.use[data.use < disp.min] <- disp.min
  
  annot.use <- intersect(colnames(object@meta.data), annot.use)
  
  if(is.list(annot.colors)){
    annot_colors <- annot.colors
  }
  if(is.null(annot.colors)){
    annot_colors <- NA
  }
  if(!is.null(pdf.file)){
    pdf(file = pdf.file, width = pdf.width, height = pdf.height, onefile = F)
  }
  pheatmap::pheatmap(data.use, cluster_cols = F, cluster_rows = F, show_colnames = F, color = colorRampPalette(colors = c("purple","black","yellow"))(11),
                     annotation_col = object@meta.data[,annot.use, drop =F], border_color = NA,
                     gaps_col = cumsum(table(object@ident)), 
                     annotation_colors = annot_colors)
  if(!is.null(pdf.file)){
    dev.off()
  }
}
bpGOEnrich <- function(pbmc.markers,org =c("mouse","human"), workers = 2){
  require(BiocParallel)
  snowparam <- SnowParam(workers = workers, type = "SOCK")
  org <- org[1]
  bplapply(as.character(sort(unique(pbmc.markers$cluster))), FUN = function(i){
    require(clusterProfiler)
    if(org == "human"){
      require(org.Hs.eg.db)
      db <- org.Hs.eg.db
    }else{
      require(org.Mm.eg.db)
      db <- org.Mm.eg.db
    }
    require(stringr)
    require(ggplot2)
    gogenes <- subset(pbmc.markers, cluster == i)$gene
    
    engo <- enrichGO(gene         = gogenes,
                     OrgDb         = db,
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
    png(filename = paste0("GO_BP results of Cluster",i,".png"), width = 600, height = 350)
    print(dotplot(engo, title = paste0("GO_BP results of Cluster",i)) + scale_y_discrete(labels = function(x) str_wrap(x, width = 60)))
    dev.off()
  }, BPPARAM = snowparam)
  
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

geneGOEnrich <- function(genes, guess.org = T, width = 600, height = 350){
  if(guess.org){
    if(sum(toupper(genes) == genes) > 0.8*length(genes))
      org <- "human"
    else
      org <- "mouse"
  }
  require(clusterProfiler)
  require(stringr)
  gogenes <- genes
  
  if(org == "human"){
    require(org.Hs.eg.db) 
    engo <- enrichGO(gene         = gogenes,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
  }
  else{
    require(org.Mm.eg.db)
    engo <- enrichGO(gene         = gogenes,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
    
  } 
  print(dotplot(engo, title = paste0("GO_BP results of Select Genes"))+ scale_y_discrete(labels = function(x) str_wrap(x, width = 60)))
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

lysisEstimate <- function(pbmc, top_n = 20, decreasing = T,add = F){
  raw_ordered <- apply(pbmc@raw.data,2,sort, decreasing = decreasing)
  raw_sum <- colSums(raw_ordered)
  raw_sum[raw_sum == 0] <- 1 # in case of NaN 
  lysis <- data.frame(ratio = colSums(raw_ordered[1:top_n, ,drop =F])/raw_sum, 
                      {tmp <- t(apply(raw_ordered[1:top_n,,drop=F], 2, range));
                       colnames(tmp) <- c("low", "high")
                       tmp
                        })
  if(add){
    pbmc@meta.data[[paste0("lysis_top",top_n,"_ratio")]] <- lysis$ratio
    pbmc@meta.data[[paste0("lysis_top",top_n,"_low")]] <- lysis$low
    pbmc@meta.data[[paste0("lysis_top",top_n,"_high")]] <- lysis$high
    return(pbmc)
  }else{
    return(lysis)
  }
}

seurat2scenic <- function(pbmc, select_ident, save.file = NULL){
  exprMat <- pbmc@data
  annot <- pbmc@meta.data[,select_ident,drop=F]
  if(is.null(save.file)){
    stop("save.file ?")
  }else{
    save(exprMat, annot, file = save.file)
  }
}

seurat2wgcna <- function(pbmc, select_ident, save.file = NULL){
  exprMat <- pbmc@data
  annot <- pbmc@meta.data[,select_ident,drop=F]
  if(is.null(save.file)){
    stop("save.file ?")
  }else{
    save(exprMat, annot, file = save.file)
  }
}

write.gsea <- function(expr, annot, filename.prefix = "gsea", add.sys.time = T){
  suffix <- ifelse(add.sys.time, Sys.time(), "")
  filename.dataset <- paste(filename.prefix, suffix, "txt",sep = ".")
  expr <- expr[,rownames(annot)]
  dataset <- cbind(rownames(expr), NA, expr)
  colnames(dataset)[1:2] <- c("Name", "Description")
  write.table(dataset, file = filename.dataset, quote = F, sep = "\t", row.names = F)
  sapply(colnames(annot), function(x){
    label <- unique(annot[,x])
    filename.pheno <- paste(filename.prefix, x, suffix, "cls",sep = ".")
    write(paste(length(annot[,x]), length(label), 1, sep = " "), file = filename.pheno)
    write(paste("#", paste(label, collapse = " "), sep = " "), file = filename.pheno, append = T)
    write(paste(annot[,x], collapse = " "), file = filename.pheno,append = T)
    #unlink(filename.pheno)
    return()
  })
  return()
}
