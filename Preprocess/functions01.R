# add `%!in%` operator
`%!in%` <- Negate(`%in%`)

# wrapper function to make volcano plots with the enhanhedvolcano package.
Make.Volcano.markers <- function(mymarkers, FCcutoff= 1, pCutoff=10e-5, mytitle=""){
  EnhancedVolcano::EnhancedVolcano(mymarkers, 
                                   x = 'avg_log2FC', 
                                   y = 'p_val_adj',
                                   titleLabSize = 10,
                                   axisLabSize = 8,
                                   lab = row.names(mymarkers),
                                   labSize = 2,
                                   subtitle = paste0("pCutOff = ", pCutoff,", log2FCcutoff = ", FCcutoff, " (cutoff FC = ", round(2^FCcutoff,1), "x)"), 
                                   subtitleLabSize = 8,
                                   cutoffLineCol = "maroon", 
                                   FCcutoff = FCcutoff, 
                                   pCutoff = pCutoff, 
                                   title = paste0("DE genes ", mytitle), 
                                   #vline = c(-0.25, 0.25), vlineType = "solid", vlineCol = "black",
                                   legendPosition = 'bottom',
                                   legendLabSize = 8,
                                   legendIconSize = 2.0,
                                   drawConnectors = TRUE,
                                   widthConnectors = 0.5
  )
}

# make the combination plots of feature and violinplot for markers specified by mymarkers
FeatureVlnPlot.markers <- function(seuratobject, use_reduction = "tsne", mymarkers, cols=NULL, group=NULL){
  intersection <- mymarkers%in%row.names(seuratobject)
  cluster.markers <- mymarkers[intersection]
  notthere <- mymarkers%!in%row.names(seuratobject)
  unavailable <- mymarkers[notthere]
  
  # make a violinplot if the intersection is not empty
    if (length(intersection)!=0){
    baseplot <- VlnPlot(seuratobject, features = cluster.markers, log = TRUE, group.by = group, cols=cols, combine = FALSE)
    baseplot <- lapply(X = baseplot, FUN = function(x) x + theme(plot.title = element_text(size = 8), axis.text=element_text(size=6),axis.title=element_text(size=6),legend.position = "none"))
    if (length(baseplot)!=2){
      k <- 0 
      for (i in (length(intersection)+1):2){#explication ici entre l'initialisation à K=0 puis une incrémentation au sein d'une boucle?
        k <- k + 1
        textplot <- ggplot() +
          theme_void() +
          annotate(geom = "text", x=1,y=1, label=paste0(unavailable[k]," is not available.")) +
          xlab(NULL)
        baseplot[[i]] <- textplot
      }
    }
  } else {
    k <- 0
    baseplot <- list()
    for (i in 1:(length(mymarkers))){#même question ici que plus haut. On remarque que k n'est jamais égal à 0 contrairement à i. Il est toujours égal à i+1
      k <- k + 1
      textplot <- eval(substitute(ggplot() +
                                    theme_void() +
                                    annotate(geom = "text", x=1,y=1, label=paste0(unavailable[k]," is not available.")) +
                                    xlab(NULL), list(i=i)))
      baseplot[[i]] <- textplot
    }
    
  }
  
  ## featureplot for the same markers
  if (length(intersection)!=0){
    fp_baseplot <- FeaturePlot(seuratobject, reduction = use_reduction, features = cluster.markers, pt.size = 1.5, label = FALSE, label.size = 2, combine = FALSE)
    fp_baseplot <- lapply(X = fp_baseplot, FUN = function(x) x + theme(plot.title = element_text(size = 8), axis.text=element_text(size=8),axis.title=element_text(size=8)))
    
    if (length(fp_baseplot)!=2){
      k <- 0
      for (i in (length(intersection)+1):2){
        k <- k + 1
        fp_textplot <- ggplot() +
          theme_void() +
          annotate(geom = "text", x=1,y=1, label=paste0(unavailable[k]," is not available.")) +
          xlab(NULL)
        fp_baseplot[[i]] <- fp_textplot
        
      }
    }
  } else {
    k <- 0
    fp_baseplot <- list()
    for (i in 1:(length(mymarkers))){
      k <- k + 1
      fp_textplot <- eval(substitute(ggplot() +
                                       theme_void() +
                                       annotate(geom = "text", x=1,y=1, label=paste0(unavailable[k]," is not available.")) +
                                       xlab(NULL), list(i=i)))
      fp_baseplot[[i]] <- fp_textplot
    }
    
  }
  
  print(wrap_plots(fp_baseplot) / wrap_plots(baseplot) + plot_annotation(title = "Marker expression per cluster"))
}

