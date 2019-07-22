# print a gene vector out as nrow * ncol
write.csv(formatGeneList(consensus_genes, nrow = 25), file = "clipboard", row.names = F, quote = F)
formatGeneList <- function(x, nrow = NULL, ncol = NULL){
  if(is.null(nrow) & is.null(ncol)) stop('specify nrow or ncol by yourself')
  if(!(is.null(nrow) | is.null(ncol))) stop('do not specify nrow and ncol concurrently')
  if(is.null(nrow)){
    nrow <- ceiling(length(x)/ncol)
  }
  if(is.null(ncol)){
    ncol <- ceiling(length(x)/nrow)
  }

  res <- list()
  for(i in 1:nrow){
    ind_start <- (i-1)*ncol +1
    ind_end <- i*ncol
    res[[as.character(i)]] <- paste0(x[ind_start:ind_end], collapse = ", ")
  }
  Reduce(f = function(x, y) {paste(x, y, sep = ",\n")}, x = res)
}

# a circle visualization for interactions between 2 cell type
circleInteractionBetween2CellType <- function(inner_r ,
                                              outer_r,
                                              pair_inner_as_source, 
                                              pair_inner_as_target, 
                                              above_horizon_degree = 5,
                                              target_r = 0.08,
                                              nodes_order = c('random','degree', "asis", "abc", 'manual'),
                                              order_random = 666,
                                              manual_order = list(up_inner = NULL,
                                                                  up_outer = NULL,
                                                                  dn_inner = NULL,
                                                                  dn_outer = NULL),
                                              curve.curvature = 0.5, curve.angle = 90, curve.ncp = 5){
  # inner_type outer_type is charactor
  # pair_inner_as_* is data.frame with 2 columns (source, target)
  # nodes_order : asis is implemented, others to be done
  # I didnot consider  genes which can be ligand and receptor stimulatory
  nodes_order <- nodes_order[1]
  if(nodes_order == "random"){
    permutateRow <- function(x){
      x[sample(1:nrow(x), nrow(x), replace = F),]
    }
    set.seed(order_random)
    pair_inner_as_target <- permutateRow(pair_inner_as_target)
    pair_inner_as_source <- permutateRow(pair_inner_as_source)
  }
  if(nodes_order == "abc"){
    pair_inner_as_target <- pair_inner_as_target[order(pair_inner_as_target[,2]),]
    pair_inner_as_source <- pair_inner_as_source[order(pair_inner_as_source[,1]),]
  }
  up_inner <- unique(pair_inner_as_target[,2])
  up_outer <- unique(pair_inner_as_target[,1])
  dn_inner <- unique(pair_inner_as_source[,1])
  dn_outer <- unique(pair_inner_as_source[,2])
  if(nodes_order == 'manual'){
    up_inner <- manual_order$up_inner
    up_outer <- manual_order$up_outer
    dn_inner <- manual_order$dn_inner
    dn_outer <- manual_order$dn_outer
  }
  if(nodes_order == 'degree'){ # needs to be optimized
    dynamicOrder <- function(net, degree, order = list(source = NULL, target = NULL)){
      if(nrow(degree) == 0){
        return(order)
      }else{
        degree_source <- subset(degree, type == 'source')
        degree_target <- subset(degree, type == 'target')
        source_list <- degree_source$gene[degree_source$degree == min(degree_source$degree)]
        target_list <- degree_target$gene[degree_target$degree == min(degree_target$degree)]
        pair_list <- unlist(net[(net[,1] %in% source_list) & (net[,2] %in% target_list),,drop = F])
        if(length(pair_list) == 0){
          if(length(source_list) == 0){
            pair_list <- target_list
          }else if(length(target_list) == 0){
            pair_list <- source_list
          }else if(length(source_list) > length(target_list)){
            source_list <- NULL
            pair_list <- target_list
          }else{
            target_list <- NULL
            pair_list <- source_list
          }
        }
        source_list <- intersect(source_list, pair_list)
        target_list <- intersect(target_list, pair_list)
        order$source <- c(order$source, source_list)
        order$target <- c(order$target, target_list)
        return(dynamicOrder(net, degree[!(degree$gene %in% c(source_list, target_list)),], order))
      }
    }
    
    getOrderByDegree <- function(net){
      source <- unique(net[,1])
      target <- unique(net[,2])
      degree <- data.frame(gene = c(source, target), 
                           type = c(rep('source',length(source)),
                                    rep('target', length(target))),
                           degree = NA)
      for(i in source){
        degree$degree[match(i, degree$gene)] <- sum(net[,1] == i)
      }
      for(i in target){
        degree$degree[match(i, degree$gene)] <- sum(net[,2] == i)
      }
      
      degree <- degree[order(degree$degree, decreasing = T),]
      return(dynamicOrder(net, degree))
    }
    res_up <- getOrderByDegree(net = pair_inner_as_target)
    res_dn <- getOrderByDegree(net = pair_inner_as_source)
    up_inner <- res_up$target
    up_outer <- res_up$source
    dn_inner <- res_dn$source
    dn_outer <- res_dn$target
  }
  
  getXY <- function(items = NULL, from = 0, to = 180, r = 1){
    deg_iv <- seq(from = from*pi/180, to = to*pi/180, length.out = length(items))
    return(data.frame(items = items, x = r*cos(deg_iv), y = r*sin(deg_iv)))
  }
  
  xy_up_inner <- getXY(up_inner, above_horizon_degree, 180 - above_horizon_degree, r = inner_r)
  xy_dn_inner <- getXY(dn_inner, above_horizon_degree + 180, 360 - above_horizon_degree, r = inner_r)
  xy_up_outer <- getXY(up_outer, above_horizon_degree, 180 - above_horizon_degree, r = outer_r)
  xy_dn_outer <- getXY(dn_outer, above_horizon_degree + 180, 360 - above_horizon_degree,  r = outer_r)
  
  getLinks <- function(source_target, xy_source, xy_target, target_r){
    source_target <- as.data.frame(source_target)
    links = as.data.frame(cbind(
      xy_source[match(source_target[,1], xy_source[,1]), c(2, 3)],
      xy_target[match(source_target[,2], xy_target[,1]), c(2, 3)]
    ))
    colnames(links) <- c("x", "y", "xend", "yend")
    
    # shrink for arrow head
    xend <- links$xend + target_r*(links$x-links$xend)/sqrt((links$y-links$yend)^2 + (links$x-links$xend)^2)
    yend <- links$yend + target_r*(links$y-links$yend)/sqrt((links$y-links$yend)^2 + (links$x-links$xend)^2)
    links$xend <- xend
    links$yend <- yend
    return(links)
  }
  
  library(ggplot2)
  library(ggrepel)
  ggplot() +
    geom_curve(mapping = aes(x = x, y = y, xend = xend, yend = yend), color = "grey",
               data = getLinks(ec_as_receptor, xy_up_outer, xy_up_inner, target_r = target_r),
               arrow = arrow(angle = 10, length = unit(0.1, "inches"),type = 'closed'),
               curvature = curve.curvature, angle = curve.angle, ncp = curve.ncp) +
    geom_point(mapping = aes(x, y), size = 8, color = 'cyan', data = xy_up_outer) +
    geom_point(mapping = aes(x, y), size = 8, color = "orange", data = xy_up_inner) +
    geom_text_repel(mapping = aes(x, y, label = items), data = xy_up_outer) +
    geom_text_repel(mapping = aes(x, y, label = items), data = xy_up_inner) +
    geom_curve(mapping = aes(x = x, y = y, xend = xend, yend = yend),  color = "grey",
               data = getLinks(ec_as_ligand, xy_dn_inner, xy_dn_outer, target_r = target_r),
               arrow = arrow(angle = 10, length = unit(0.1, "inches"),type = 'closed'),
               curvature = curve.curvature, angle = curve.angle, ncp = curve.ncp) +
    geom_point(mapping = aes(x, y), size = 8, color = "orange", data = xy_dn_inner) +
    geom_point(mapping = aes(x, y), size = 8, color = 'cyan', data = xy_dn_outer) +
    geom_text_repel(mapping = aes(x, y, label = items), data = xy_dn_inner) +
    geom_text_repel(mapping = aes(x, y, label = items), data = xy_dn_outer) +
    geom_hline(yintercept = 0) +
#    coord_fixed(ratio = 1) +
    theme_void()
}

circleInteractionBetween2CellType(inner_r = 1.3, outer_r = 2,  above_horizon_degree = 10,
                                  pair_inner_as_source = as.data.frame(ec_as_ligand), 
                                  pair_inner_as_target = as.data.frame(ec_as_receptor),
                                  curve.curvature = 0.3, curve.angle = 30, curve.ncp = 20,
                                  nodes_order = 'random',  order_random = 666)
