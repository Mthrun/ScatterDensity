DensityScatter.DDCAL = function (X,
                                 Y,
                                 nClusters = 12,
                                 Plotter = "native",
                                 SDHorPDE = TRUE,
                                 PDEsample = 10000,
                                 Marginals = FALSE,
                                 na.rm = TRUE,
                                 pch = 10,
                                 Size = 1,
                                 xlab="x",
                                 ylab="y",
                                 main = "",
                                 lwd = 2,
                                 xlim=NULL,
                                 ylim=NULL,
                                 Polygon,
                                 BW = TRUE,
                                 Silent = FALSE,
                                 ...) {
  # 
  if(length(X) != length(Y)){
    stop("DensityScatter.DDCAL: length of X points does not equal lengt of Y-points")
  }
  if(missing(xlab)){
    xlab = deparse1(substitute(X))
  }
  if(missing(ylab)){
    ylab = deparse1(substitute(Y))
  }
  if(isTRUE(na.rm)){
    ind.na = is.finite(X) & is.finite(Y)
    X      = X[ind.na]
    Y      = Y[ind.na]
  }
  Bool = FALSE
  if(length(X) > PDEsample && isFALSE(SDHorPDE)) {
    Bool = TRUE
    ind  = sample(1:length(X), PDEsample)
  }else{
    ind  = 1:length(X)
  }
  if(isFALSE(Silent)){
    message("DensityScatter.DDCAL: Estimating Density...")
  }
  if(isTRUE(SDHorPDE)){
    Dens = SmoothedDensitiesXY(X = X[ind], Y = Y[ind], PlotIt = FALSE)$Densities
  }else{
    Dens = PDEscatter(x = X[ind], y = Y[ind], PlotIt = -1)$Densities
  }
  
  if(isFALSE(Silent)){
    message("DensityScatter.DDCAL: Estimating colors...")
  }
  
  cls = DDCAL(Dens, nClusters = nClusters, minBoundary = 0.2, 
              maxBoundary = 0.6, numSimulations = 20, 
              csTolerance = 0.45, csToleranceIncrease = 0.5)

	if(missing(xlim)){
    xlim = c(min(X, na.rm = T), max(X, na.rm = T))
	}
	
	if(missing(ylim)){
		ylim = c(min(Y, na.rm = T), max(Y, na.rm = T))
	}
  if(!missing(Polygon)){
    Polygon=assertPolygon(Polygon)
  }
  
  if(Plotter %in% c("native", "ggplot2")){
    cols = rep(NA_character_, length(ind))
    ncolors  = 1
    ncolors2 = 1
    colpalette = colorRampPalette(c("navyblue", "darkblue", rep("blue", 
                                                                ncolors2), rep("turquoise", ncolors2), rep("green", ncolors), 
                                    rep("chartreuse", ncolors), rep("yellow", ncolors), rep("orange", 
                                                                                            ncolors), rep("red", ncolors), rep("darkred", ncolors2)))
    
    colpal <- cut(cls, length(cls), labels = FALSE)
    cols <- colpalette(length(cls))[colpal]
  }
  
  if (Plotter == "native") {
    if (isTRUE(Marginals)) {
      def.par <- par(no.readonly = TRUE)
      on.exit(par(def.par))
      
      m <- graphics::layout(matrix(c(1, 2, 2, 2, 1, 2, 
                                     2, 2, 1, 2, 2, 2, 4, 3, 3, 3), 4, 4))
      par(fig = c(0, 0.8, 0, 0.8), 
          mar = c(5, 4, 0, 0))
      # If marginal plots are used main should be plotted above the x-Axis marginal plot instead the scatter plot
      mainMarginal = main
      main = ""
    }
    if (isFALSE(Silent)) 
      message("DensityScatter.DDCAL: Plotting...")
    if (isTRUE(Bool)) {
      if (length(X) > 10 * PDEsample) {
        ind2 = sample(1:length(X), 3 * PDEsample)
      }
      else {
        ind2 = 1:length(X)
      }
      plot(X[ind2], Y[ind2], xlab = xlab, ylab = ylab, 
           xlim = xlim, ylim = ylim, col = "navyblue", pch = 3, 
           cex = Size,main=main, ...)
      points(x = X[ind], y = Y[ind], col = cols, pch = pch, 
             cex = Size, ...)
    }
    else {
      plot(x = X, y = Y, xlab = xlab, ylab = ylab, xlim = xlim, 
           ylim = ylim, col = cols, pch = pch, cex = Size,main=main, ...)
    }
    if(!missing(Polygon)){
      points(Polygon[,1],Polygon[,2],lwd=3,col="magenta",type="l")
    }
    
    if (isTRUE(Marginals)) {
      par(fig = c(0, 0.8, 0.8, 
                  1), mar = c(0, 4, 4, 0), new = TRUE)
      V = DataVisualizations::ParetoDensityEstimation(X[ind], PlotIt = F)
      plot(V$kernels, V$paretoDensity, xlim = xlim, ylim = c(0, max(V$paretoDensity)), 
           ylab = "PDE", type = "l", lwd = lwd, axes = FALSE, bty = "n", main = mainMarginal)
      axis(side = 2, at = pretty(V$paretoDensity, n = 3), 
           labels = gsub("\\.?0*$", "", formatC(pretty(V$paretoDensity, n = 3), format = "f")))
      par(fig = c(0.8, 1, 0, 0.8), 
          mar = c(5, 0, 0, 2), new = TRUE)
      V = DataVisualizations::ParetoDensityEstimation(Y[ind], PlotIt = F)
      plot(V$paretoDensity, V$kernels, xlim = c(0, max(V$paretoDensity)), 
           ylim = ylim, xlab = "PDE", type = "l", lwd = lwd, axes = FALSE, bty = "n")
      axis(side = 1, at = pretty(V$paretoDensity, n = 3), 
           labels = gsub("\\.?0*$", "", formatC(pretty(V$paretoDensity, n = 3), format = "f")))
      
      #par(def.par)#apparantly on.exit is the better solution
    }
  }else if (Plotter == "ggplot2") {
    if (isFALSE(Silent)) 
      message("DensityScatter.DDCAL: Preparing for ggplot2..")
    if (isTRUE(Bool)) {
      if (length(X) > 2 * PDEsample) {
        mod = ClusterR::MiniBatchKmeans(cbind(X, Y), 
                                        2 * PDEsample)
        Centroids = mod$centroids
        DFnull = data.frame(x2 = Centroids[, 1], y2 = Centroids[,2], Colors = "navyblue")
      }
      else {
        ind2 = 1:length(X)
        DFnull = data.frame(x2 = X, y2 = Y, Colors = "navyblue")[ind2,]
      }
      if (isFALSE(Silent)) 
        message("DensityScatter.DDCAL: using ggplot2..")
      DF = data.frame(x = X[ind], y = Y[ind], Colors = cols, 
                      Cls = cls)
      ggobj = ggplot2::ggplot(DF, ggplot2::aes_string(x = "x", y = "y", 
                                               col = "Colors")) + ggplot2::geom_point(data = DFnull, 
                                                                                      mapping = ggplot2::aes_string(x = "x2", y = "y2", col = "Colors"), 
                                                                                      size = Size, shape = 3, alpha = 0.4) + 
        ggplot2::geom_point(size = Size, shape = pch, 
                            alpha = 0.4) + ggplot2::theme(legend.position = "none") + 
        ggplot2::scale_color_identity()
    }
    else {
      DF = data.frame(x = X, y = Y, Colors = cols, Cls = cls)
      ggobj = ggplot2::ggplot(DF, ggplot2::aes_string(x = "x", y = "y", 
                                               col = "Colors"), alpha = 0.05) +
        ggplot2::geom_point(size = Size, shape = pch) + 
        ggplot2::theme(legend.position = "none") + ggplot2::scale_color_identity()
    }
    if(!missing(Polygon)) {
      DFPoly = data.frame(PolygonX = Polygon[,1], PolygonY = Polygon[,2])
      ggobj = ggobj + ggplot2::geom_polygon(data = DFPoly, ggplot2::aes_string(x = "PolygonX", y = "PolygonY"), 
                                            fill = "transparent", color = "magenta", size = 2)
    }
    ggobj = ggobj + xlab(xlab) + ylab(ylab)
    if (isTRUE(BW)) 
      ggobj = ggobj + ggplot2::theme_bw()
    if(nchar(main) != 0)
      ggobj = ggobj + ggplot2::ggtitle(main) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    if (isTRUE(Marginals)) {
      ggobj2 <- ggExtra::ggMarginal(ggobj, type = "density")
      #print(ggobj2)
      return(ggobj2)
    }
    else {
      #print(ggobj)
      return(ggobj)
    }
  }else{
    if(isFALSE(Silent)){
      message("DensityScatter.DDCAL: Preparing plotly visualization...")
    }
    
    palette = colorRampPalette(c("darkblue", "blue", "lightblue1", "green","yellow", "red", "darkred"))
    
    cls = cls + 1
    Colors  = palette(length(unique(cls)))
    
    if(is.null(xlim)){
      xlim = c(min(X), max(X)+1)
    }
    if(is.null(ylim)){
      ylim = c(min(Y), max(Y)+2)
    }
    
    if(isFALSE(Marginals)){
      
      #mod = ClusterR::MiniBatchKmeans(cbind(X, Y), 2 * PDEsample)
      #x = Sys.time()
      #mod = ClusterR::MiniBatchKmeans(data = cbind(X, Y), clusters = 2 * PDEsample)
      #VKM = FCPS::kmeansClustering(cbind(X,Y), ClusterNo = 2*PDEsample)
      #y = Sys.time()
      #y-x
      # BatchKmeans
      # 10k samples =>  2.124521 mins
      
      # 10k Samples => 4.951251 mins
      #  5k Samples => 5.877196 mins
      #  2k Samples => 1.521917 mins
      
      #Centroids = mod$centroids
      #Idx = sample(1:length(X), 10000, replace = FALSE)
      #KMeansCls = VKM$Cls
      #CentroidsIdx = unlist(lapply(unique(KMeansCls), function(x, KMeansCls){
      #  which(KMeansCls == x)[1]
      #}, KMeansCls))
      #Idx = as.numeric(CentroidsIdx)
      
      plotOut = plotly::plot_ly()
      #plotOut = plotly::add_markers(p = plotOut,
      #                              x = X[Idx], y = Y[Idx],
      #                              marker = list(color = Colors[cls[Idx]], size = 3), type = "scatter")
      plotOut = plotly::add_markers(p = plotOut,
                                    x = X[ind], y = Y[ind],
                                    marker = list(color = Colors[cls], size = 3), type = "scatter")
      if(!missing(Polygon)){
        plotOut = plotly::add_polygons(p = plotOut,
                                       x = Polygon[,1], y = Polygon[,2], color = I("magenta"), opacity = 0.7)
      }
      #plotOut = plotly::add_lines(p = plotOut,
      #                            x = Polygon[,1], y = Polygon[,2],
      #                            marker = list(color = Colors[cls], size = 3), type = "scatter")
      
      plotOut = plotly::layout(p      = plotOut,
                               title  = main,
                               xaxis  = list(title = xlab, fixedrange = T),#, scaleanchor="y", scaleratio=1),
                               yaxis  = list(title = ylab, fixedrange = T),
                               plot_bgcolor = "rgb(254, 254, 254)",              # plot_bgcolor = "rgb(254, 247, 234)",
                               paper_bgcolor = "rgb(254, 254, 254)")             # paper_bgcolor = "rgb(254, 247, 234)"
      plotOut = plotly::hide_colorbar(p = plotOut)
      plotOut = plotly::hide_legend(p = plotOut)
      plotOut = plotly::config(p = plotOut, displayModeBar=T, editable=T)
      print(plotOut)
    }else{
      PDE1 = DataVisualizations::PDEplot(Data = X[ind])
      PDE2 = DataVisualizations::PDEplot(Data = Y[ind])
      DomX = PDE1$kernels
      ValX = PDE1$paretoDensity
      DomY = PDE2$kernels
      ValY = PDE2$paretoDensity
      
      Fig1 = plotly::plot_ly(x = DomX, y = ValX, type = 'scatter', mode = "lines",
                     alpha =.5, line = list(color = "black"))
      Fig2 = plotly::plotly_empty()
      Fig3 = plotly::plot_ly(x = X[ind], y = Y[ind], type = 'scatter',
                             mode = 'markers', alpha = .5,
                             marker = list(color = Colors[cls], size = 3))
      if(!missing(Polygon)){
        Fig3 = plotly::add_polygons(p = Fig3,
                                    x = Polygon[,1], y = Polygon[,2], color = I("magenta"), opacity = 0.7)
      }
      #layout(p = Fig1, yaxis = list(showline = TRUE), xaxis = list(showline = TRUE))
      Fig4 = plotly::plot_ly(y = DomY, x = ValY, type = 'scatter', mode = "lines",
                             alpha = .5, line = list(color = "black"))
      plotOut = plotly::subplot(Fig1, Fig2, Fig3, Fig4, nrows = 2,
                                heights = c(.2, .8), widths = c(.8,.2),
                                margin = 0, shareX = TRUE, shareY = TRUE,
                                titleX = FALSE, titleY = FALSE)
      plotOut = plotly::layout(p = plotOut, title = main, showlegend = FALSE, barmode = 'overlay',
                               yaxis = list(title = ylab, showline = FALSE),
                               xaxis = list(title = xlab, showline = FALSE))
      print(plotOut)
    }
    return(list("Plot" = plotOut, "Cls" = cls))
    #if(isFALSE(Silent)) 
    #  message("DensityScatter.DDCAL: plotly not implemented")
  }
}
