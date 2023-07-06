## SOurces :
## Sqjin. (n.d.). Arrow missing from netVisual_chord_gene · Issue #266 · sqjin/CellChat. GitHub. 
## Retrieved June 30, 2023, from https://github.com/sqjin/CellChat/issues/266
## 
## author : pedriniedoardo

library(circlize)

netVisual_chord_gene2 <- function (object, slot.name = "net", color.use = NULL, signaling = NULL, 
                                   pairLR.use = NULL, net = NULL, sources.use = NULL, targets.use = NULL, 
                                   lab.cex = 0.8, small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03), 
                                   link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, 
                                   reduce = -1, transparency = 0.4, link.border = NA, title.name = NULL, 
                                   legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE, 
                                   thresh = 0.05, ...) 
{
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }
  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }
  if (is.null(net)) {
    prob <- slot(object, "net")$prob
    pval <- slot(object, "net")$pval
    prob[pval > thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source", "target", "interaction_name")
    pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", 
                                              "pathway_name", "ligand", "receptor", "annotation", 
                                              "evidence"))
    idx <- match(net$interaction_name, rownames(pairLR))
    temp <- pairLR[idx, ]
    net <- cbind(net, temp)
  }
  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], 
                                 pairLR.use = object@LR$LRsig, key = "pathway_name", 
                                 matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }
  if (!is.null(pairLR.use)) {
    if ("interaction_name" %in% colnames(pairLR.use)) {
      net <- subset(net, interaction_name %in% pairLR.use$interaction_name)
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
    }
  }
  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source", "target", "pathway_name", 
                                "prob"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    net <- net %>% dplyr::group_by(source_target, pathway_name) %>% 
      dplyr::summarize(prob = sum(prob))
    a <- stringr::str_split(net$source_target, "sourceTotarget", 
                            simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net$ligand <- net$pathway_name
    net$receptor <- " "
  }
  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) {
      sources.use <- levels(object@idents)[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  else {
    sources.use <- levels(object@idents)
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) {
      targets.use <- levels(object@idents)[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  else {
    targets.use <- levels(object@idents)
  }
  df <- subset(net, prob > 0)
  if (nrow(df) == 0) {
    stop("No signaling links are inferred! ")
  }
  if (length(unique(net$ligand)) == 1) {
    message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
  }
  df$id <- 1:nrow(df)
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(" ", 
                                                       j - 1), collapse = ""))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(" ", 
                                                           j - 1), collapse = ""))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }
  cell.order.sources <- levels(object@idents)[levels(object@idents) %in% 
                                                sources.use]
  cell.order.targets <- levels(object@idents)[levels(object@idents) %in% 
                                                targets.use]
  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]
  order.source <- unique(df.ordered.source[, c("ligand", "source")])
  order.target <- unique(df.ordered.target[, c("receptor", 
                                               "target")])
  order.sector <- c(order.source$ligand, order.target$receptor)
  if (is.null(color.use)) {
    color.use = scPalette(nlevels(object@idents))
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source, 
                                                                         df$target))]
  }
  else if (is.null(names(color.use))) {
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source, 
                                                                         df$target))]
  }
  edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)
  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector
  df.plot <- df.ordered.source[, c("ligand", "receptor", "prob")]
  if (directional == 2) {
    link.arr.type = "triangle"
  }
  else {
    link.arr.type = "big.arrow"
  }
  circos.clear()
  chordDiagram(df.plot, order = order.sector, col = edge.color, 
               grid.col = grid.col, transparency = transparency, link.border = link.border, 
               directional = directional, direction.type = c("diffHeight+arrows"), link.arr.type = link.arr.type, annotationTrack = "grid", 
               annotationTrackHeight = annotationTrackHeight, preAllocateTracks = list(track.height = max(strwidth(order.sector))), 
               small.gap = small.gap, big.gap = big.gap, link.visible = link.visible, 
               scale = scale, link.target.prop = link.target.prop, reduce = reduce, 
               ...)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                niceFacing = TRUE, adj = c(0, 0.5), cex = lab.cex)
  }, bg.border = NA)
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(color.use), 
                                  type = "grid", legend_gp = grid::gpar(fill = color.use), 
                                  title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc") - unit(legend.pos.x, 
                                                        "mm"), y = unit(legend.pos.y, "mm"), just = c("right", 
                                                                                                      "bottom"))
  }
  circos.clear()
  if (!is.null(title.name)) {
    text(-0, 1.02, title.name, cex = 1)
  }
  gg <- recordPlot()
  return(gg)
}