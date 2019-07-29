options(stringsAsFactors = F)
rm(list = ls())

# given source-target network 
# output order of source and target for layout
getForceOrder <- function(net, iter = 100, seed = 666){
  source <- unique(net[,1])
  target <- unique(net[,2])
  nodes <- c(source, target)
  set.seed(seed)
  xy <- data.frame(x = sample(x = 1:length(nodes), size = length(nodes), replace = F),
                   y = 1, row.names = nodes)
  xy[target, 'y'] <- -1 # target_y = -1, source_y = 1
  
  kk <- (max(xy[,1]) - min(xy[,1]))/length(nodes)
  ka <- 100*kk
  kr <- kk
  max_disp <- 8
  # fa <- function(d, k){return(d^2/k)}
  # fr <- function(d, k){return(k^2/d)}
  mod <- function(d){abs(d[1])}
  for(i in 1:iter){
    # initial pos displacement
    u_disp = data.frame(x = rep(0, length(nodes)),
                        y = rep(0, length(nodes)), 
                        row.names = nodes)
    # repulsive 
    for(i in source){
      for(j in source){
        if(i == j) next()
        d <- xy[i,1] - xy[j,1]
        if(d == 0) d <- 1 # in case d==0
        u_disp[i,1] = u_disp[i,1] + (kr^2)/d
      }
    }
    for(i in target){
      for(j in target){
        if(i == j) next()
        d <- xy[i,1] - xy[j,1]
        if(d == 0) d <- 1 # in case d==0
        u_disp[i,1] = u_disp[i,1] + (kr^2)/d
      }
    }
    # attraction
    for(k in 1:nrow(net)){
      u <- net[k,1]
      v <- net[k,2]
      d <- xy[v,1] - xy[u,1]
      u_disp[u,1] = u_disp[u,1] + sign(d)*(d^2)/ka
      u_disp[v,1] = u_disp[v,1] - sign(d)*(d^2)/ka
    }
    
    for(l in nodes){
      u_disp_l <- mod(u_disp[l,])
      if(u_disp_l > max_disp)
        u_disp[l,1] <- u_disp[l,1]*max_disp/u_disp_l
      xy[l,1] = xy[l,1] + u_disp[l,1]
    }
  }
  #return(xy)
  return(list(source[order(xy[source,1])], 
              target[order(xy[target,1])]))
}

# Plot source-target network as the top and bottom semocircles
circleInteractionBetween2CellType <- function(pair_inner_as_source, 
                                              pair_inner_as_target, 
                                              above_horizon_degree = c(5, 10),
                                              inner_r ,
                                              outer_r,
                                              size = 1, target_r = 0.01*size,
                                              nodes_order = c('force-directed','manual', 'random', "asis", "abc"),
                                              random_seed = 666,
                                              force_iter = 100,
                                              manual_order = list(up_inner = NULL,
                                                                  up_outer = NULL,
                                                                  dn_inner = NULL,
                                                                  dn_outer = NULL),
                                              link.curvature = 0.5, link.angle = 90, link.ncp = 5){
  # inner_type outer_type is charactor
  # pair_inner_as_* is data.frame with 2 columns (source, target)
  # nodes_order : asis is implemented, others to be done
  # I didnot consider  genes which can be ligand and receptor stimulatory
  nodes_order <- nodes_order[1]
  if(nodes_order == "random"){
    permutateRow <- function(x){
      x[sample(1:nrow(x), nrow(x), replace = F),]
    }
    set.seed(random_seed)
    pair_inner_as_target <- permutateRow(pair_inner_as_target)
    pair_inner_as_source <- permutateRow(pair_inner_as_source)
  }
  if(nodes_order == "abc"){
    pair_inner_as_target <- pair_inner_as_target[order(pair_inner_as_target[,2], pair_inner_as_target[,1]),]
    pair_inner_as_source <- pair_inner_as_source[order(pair_inner_as_source[,1], pair_inner_as_source[,2]),]
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
  
  if(nodes_order == 'force_directed'){
    tmp <- getForceOrder(net = pair_inner_as_target, iter = force_iter, seed = random_seed)
    up_inner <- tmp[[2]]
    up_outer <- tmp[[1]]
    tmp <- getForceOrder(net = pair_inner_as_source, iter = force_iter, seed = random_seed)
    dn_inner <- tmp[[1]]
    dn_outer <- tmp[[2]]
  }
  
  getXY <- function(items = NULL, from = 0, to = 180, r = 1){
    deg_iv <- seq(from = from*pi/180, to = to*pi/180, length.out = length(items))
    return(data.frame(items = items, x = r*cos(deg_iv), y = r*sin(deg_iv)))
  }
  
  xy_up_inner <- getXY(up_inner, above_horizon_degree[1], 180 - above_horizon_degree[1], r = inner_r)
  xy_dn_inner <- getXY(dn_inner, above_horizon_degree[1] + 180, 360 - above_horizon_degree[1], r = inner_r)
  xy_up_outer <- getXY(up_outer, above_horizon_degree[2], 180 - above_horizon_degree[2], r = outer_r)
  xy_dn_outer <- getXY(dn_outer, above_horizon_degree[2] + 180, 360 - above_horizon_degree[2],  r = outer_r)
  
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
               data = getLinks(pair_inner_as_target, xy_up_outer, xy_up_inner, target_r = target_r),
               arrow = arrow(angle = 20, length = unit(0.1, "inches"),type = 'closed'),
               curvature = -link.curvature, angle = link.angle, ncp = link.ncp) +
    geom_point(mapping = aes(x, y), size = size, color = 'cyan', data = xy_up_outer) +
    geom_point(mapping = aes(x, y), size = size, color = "orange", data = xy_up_inner) +
    geom_text_repel(mapping = aes(x, y, label = items), data = xy_up_outer) +
    geom_text_repel(mapping = aes(x, y, label = items), data = xy_up_inner) +
    geom_curve(mapping = aes(x = x, y = y, xend = xend, yend = yend),  color = "grey",
               data = getLinks(pair_inner_as_source, xy_dn_inner, xy_dn_outer, target_r = target_r),
               arrow = arrow(angle = 20, length = unit(0.1, "inches"),type = 'closed'),
               curvature = -link.curvature, angle = link.angle, ncp = link.ncp) +
    geom_point(mapping = aes(x, y), size = size, color = "orange", data = xy_dn_inner) +
    geom_point(mapping = aes(x, y), size = size, color = 'cyan', data = xy_dn_outer) +
    geom_text_repel(mapping = aes(x, y, label = items), data = xy_dn_inner) +
    geom_text_repel(mapping = aes(x, y, label = items), data = xy_dn_outer) +
    geom_hline(yintercept = 0) +
    #    coord_fixed(ratio = 1) +
    theme_void()
}

# generage random source-target network for illustration
generate_net <- function(source, target, n_links = round(sqrt(length(source) * length(target))*3), seed = 666){
  set.seed(seed)
  return(unique(data.frame(source = sample(source, size = n_links, replace = T),
                           target = sample(target, size = n_links, replace = T))))
}
up_net <- generate_net(paste0("MC-S", 1:18), paste0("EC-T", 1:22), seed = 111)
dn_net <- generate_net(paste0("EC-S", 1:18), paste0("MC-T", 1:22), seed = 222)

circleInteractionBetween2CellType(inner_r = 1.2, outer_r = 2,  above_horizon_degree = c(10,10), 
                                  pair_inner_as_source = up_net, 
                                  pair_inner_as_target = dn_net,
                                  link.curvature = -0.1, link.angle = 30, link.ncp = 20,size = 6,
                                  nodes_order = 'random',  random_seed = 888)
circleInteractionBetween2CellType(inner_r = 1.2, outer_r = 2,  above_horizon_degree = c(10,10), 
                                  pair_inner_as_source = up_net, 
                                  pair_inner_as_target = dn_net,
                                  link.curvature = -0.1, link.angle = 30, link.ncp = 20,size = 6,
                                  nodes_order = 'force_directed',  
                                  force_iter = 100,
                                  random_seed = 888)
