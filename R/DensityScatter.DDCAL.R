DensityScatter.DDCAL = function (X,
                                 Y,
                                 nClusters = 12,
                                 Plotter = "native",
                                 SDHorPDE = TRUE,
                                 LimitShownPoints=FALSE,
                                 Marginals = FALSE,
                                 na.rm = TRUE,
                                 pch,
                                 Size,
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
  
  if(isFALSE(Silent)){
    message("DensityScatter.DDCAL: Estimating Density...")
  }
  if(isTRUE(SDHorPDE)){
    Dens = SmoothedDensitiesXY(X = X, Y = Y, PlotIt = FALSE,Compute="Parallel")$Densities
  }else{
    Dens = PDEscatter(x = X, y = Y, PlotIt = -1,Compute="Parallel")$Densities
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
  
    cols = rep(NA_character_, length(X))
    ncolors  = 1
    ncolors2 = 1
    colpalette = colorRampPalette(c("navyblue", "darkblue", rep("blue", 
                                                                ncolors2), rep("turquoise", ncolors2), rep("green", ncolors), 
                                    rep("chartreuse", ncolors), rep("yellow", ncolors), rep("orange", 
                                                                                            ncolors), rep("red", ncolors), rep("darkred", ncolors2)))
    
    colpal <- cut(cls, length(cls), labels = FALSE)
    cols <- colpalette(length(cls))[colpal]
    
    ##Set defaulgt plotting parameters
    #nota: if I do first cbind then as.data.frame() the x and y ticks are overburden with numbers in plotly and ggplot2
    #dont know why
    switch (Plotter,
            "native" = {
              DF = cbind(x = X, y = Y, Colors = cols, Cls = cls)
              if(missing(Size)) Size=1
              if(missing(pch)) pch=20#10
              
              if (isFALSE(Marginals)) {
                bty="n"
              }else{
                bty="o"
              }
            },
            "ggplot2" = {
              DF = data.frame(x = X, y = Y, Colors = cols, Cls = cls)
              opacity = 0.05 #alpha
              if(missing(Size)) Size=3
              if(missing(pch)) pch=20#10
            },
            "plotly" = {
              DF = data.frame(x = X, y = Y, Colors = cols, Cls = cls)
              opacity = 0.85
              if(missing(Size)) Size=6
              if(missing(pch)) pch=0#127
            },
            {
              DF = cbind(x = X, y = Y, Colors = cols, Cls = cls)
            }
    )
    
    if (isTRUE(LimitShownPoints)) {
      subsample_ind=SampleScatter(X,Y,ThresholdPoints = 10)
      DF=DF[subsample_ind,]
    }

  ##Plotting ----
  ## native plot ----
    switch (Plotter,
            "native" = {
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

      plot(x = DF[,"x"], y = DF[,"y"], xlab = xlab, ylab = ylab, xlim = xlim, 
           ylim = ylim, col = DF[,"Colors"], pch = pch, cex = Size,main=main,bty = bty, ...)
    
    
    if(!missing(Polygon)){
      points(Polygon[,1],Polygon[,2],lwd=3,col="magenta",type="l")
    }
    
    if (isTRUE(Marginals)) {
      par(fig = c(0, 0.8, 0.8, 
                  1), mar = c(0, 4, 4, 0), new = TRUE)
      V = DataVisualizations::ParetoDensityEstimation(X, PlotIt = F)
      plot(V$kernels, V$paretoDensity, xlim = xlim, ylim = c(0, max(V$paretoDensity)), 
           ylab = "PDE", type = "l", lwd = lwd, axes = FALSE, bty = "n", main = mainMarginal)
      axis(side = 2, at = pretty(V$paretoDensity, n = 3), 
           labels = gsub("\\.?0*$", "", formatC(pretty(V$paretoDensity, n = 3), format = "f")))
      par(fig = c(0.8, 1, 0, 0.8), 
          mar = c(5, 0, 0, 2), new = TRUE)
      V = DataVisualizations::ParetoDensityEstimation(Y, PlotIt = F)
      plot(V$paretoDensity, V$kernels, xlim = c(0, max(V$paretoDensity)), 
           ylim = ylim, xlab = "PDE", type = "l", lwd = lwd, axes = FALSE, bty = "n")
      axis(side = 1, at = pretty(V$paretoDensity, n = 3), 
           labels = gsub("\\.?0*$", "", formatC(pretty(V$paretoDensity, n = 3), format = "f")))
      
      #par(def.par)#apparantly on.exit is the better solution
    }
      plotOut=NULL
  }, ##ggplot2 ----
  "ggplot2"= {
    if (isFALSE(Silent)) 
      message("DensityScatter.DDCAL: Preparing for ggplot2..")
      # DF=as.data.frame(DF)
      # rownames(DF)=NULL
      ggobj = ggplot2::ggplot(DF, ggplot2::aes_string(x = "x", y = "y", 
                                                      col = "Colors"), alpha = opacity) +
        ggplot2::geom_point(size = Size, shape = pch) + 
        ggplot2::theme(legend.position = "none") + ggplot2::scale_color_identity()

    if(!missing(Polygon)) {
      DFPoly = data.frame(PolygonX = Polygon[,1], PolygonY = Polygon[,2])
      ggobj = ggobj + ggplot2::geom_polygon(data = DFPoly, ggplot2::aes_string(x = "PolygonX", y = "PolygonY"), 
                                            fill = "transparent", color = "magenta", size = 2)
    }
    
    ggobj = ggobj + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    
    if (isTRUE(BW)) {
      ggobj = ggobj + ggplot2::theme_bw()
    }
    if(nchar(main) != 0){
      ggobj = ggobj + ggplot2::ggtitle(main) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    }
    
    if (isTRUE(Marginals)) {
      ggobj <- ggExtra::ggMarginal(ggobj, type = "density",margins = "both")
    }else{
      ggobj=ggobj+ggplot2::theme(panel.border = element_blank())
    }
      print(ggobj)
      plotOut=ggobj
  },##plotly plot ----
  "plotly"= {
    if(isFALSE(Silent)){
      message("DensityScatter.DDCAL: Preparing plotly visualization...")
    }

    if(is.null(xlim)){
      xlim = c(min(X), max(X)+1)
    }
    if(is.null(ylim)){
      ylim = c(min(Y), max(Y)+2)
    }
    # DF=as.data.frame(DF)
    # rownames(DF)=NULL
    #generate main plot
    Fig3=plotly::plot_ly(
      data = DF,
      x = ~x,
      y = ~y,
      type = "scatter",
      mode = "markers",
      # Directly pass the color column, and set alpha (opacity).
      marker = list(
        color = ~Colors,
        size = Size,     # from your DF or a numeric value
        opacity = opacity,
        symbol = pch,  # Use this if you want to map shapes,
        line = list(width = 0)
        # but must remap R's pch to Plotly's symbol codes
      )
    )
    if(isFALSE(Marginals)){
      #plotOut consists of only one plot
      plotOut=Fig3
      if(!missing(Polygon)){
        plotOut = plotly::add_polygons(p = plotOut,
                                       x = Polygon[,1], y = Polygon[,2], color = I("magenta"), opacity = 0.7)
      }

      plotOut = plotly::layout(p      = plotOut,
                               title  = main,
                               xaxis  = list(title = xlab, fixedrange = T,      
                                             showline = FALSE,    # Remove x-axis line
                                             zeroline = FALSE#,    # Remove zero line
                                             #showgrid = FALSE     # Remove grid lines
                                             ),#, scaleanchor="y", scaleratio=1),
                               yaxis  = list(title = ylab, fixedrange = T,
                                             showline = FALSE,    # Remove x-axis line
                                             zeroline = FALSE    # Remove zero line
                                             #showgrid = FALSE     # Remove grid lines
                                             ),
                               plot_bgcolor = "rgb(254, 254, 254)",              # plot_bgcolor = "rgb(254, 247, 234)",
                               paper_bgcolor = "rgb(254, 254, 254)")             # paper_bgcolor = "rgb(254, 247, 234)"
      plotOut = plotly::hide_colorbar(p = plotOut)
      plotOut = plotly::hide_legend(p = plotOut)
      plotOut = plotly::config(p = plotOut, displayModeBar=T, editable=T)
      print(plotOut)
    }else{#Marginals=TRUE
      #plotOut consists of four plots of which one is empty, third one is generated above
      PDE1 = DataVisualizations::PDEplot(Data = X)
      PDE2 = DataVisualizations::PDEplot(Data = Y)
      DomX = PDE1$kernels
      ValX = PDE1$paretoDensity
      DomY = PDE2$kernels
      ValY = PDE2$paretoDensity
      
      Fig1 = plotly::plot_ly(x = DomX, y = ValX, type = 'scatter', mode = "lines",
                     alpha =.5, line = list(color = "black"))
      Fig2 = plotly::plotly_empty(type = "scatter",mode = "markers")

      if(!missing(Polygon)){
        Fig3 = plotly::add_polygons(p = Fig3,
                                    x = Polygon[,1], y = Polygon[,2], color = I("magenta"), opacity = 0.7)
      }

      Fig4 = plotly::plot_ly(y = DomY, x = ValY, type = 'scatter', mode = "lines",
                             alpha = .5, line = list(color = "black"))
      
      Fig1 = plotly::layout(p = Fig1, showlegend = FALSE,
                               yaxis = list(title = FALSE, showline = T,zeroline = FALSE),
                               xaxis = list(title = FALSE, showline = T,zeroline = FALSE))
      Fig4 = plotly::layout(p = Fig4, showlegend = FALSE,
                            yaxis = list(title = FALSE, showline = T,zeroline = FALSE),
                            xaxis = list(title = FALSE, showline = T,zeroline = FALSE))
      
      plotOut = plotly::subplot(Fig1, Fig2, Fig3, Fig4, nrows = 2,
                                heights = c(.2, .8), widths = c(.8,.2),
                                margin = 0, shareX = T, shareY = T,
                                titleX = F, titleY = F)
      plotOut = plotly::layout(p = plotOut, title = main, showlegend = FALSE, barmode = 'overlay',
                               yaxis = list(title = ylab, showline = T,zeroline = FALSE),
                               xaxis = list(title = xlab, showline = T,zeroline = FALSE))
      print(plotOut)
    }

  },{
    warning("DensityScatter.DDCAL: Returning Density as no correct plotter selected")
    plotOut=NULL
  }
  )
    
    if (isTRUE(LimitShownPoints)) {
      Dens=Dens[subsample_ind]
    }
    DF=cbind("Density" = Dens,DF)
    DF=as.matrix(DF)
    
  return(invisible(list("DF" = as.matrix(DF), "PlotHandle" = plotOut)))
}


##HALDE ----

# if (length(X) > 2 * PDEsample) {
#   mod = ClusterR::MiniBatchKmeans(cbind(X, Y), 
#                                   2 * PDEsample)
#   Centroids = mod$centroids
#   DFnull = data.frame(x2 = Centroids[, 1], y2 = Centroids[,2], Colors = "navyblue")
# }else {
#   ind2 = 1:length(X)
#   DFnull = data.frame(x2 = X, y2 = Y, Colors = "navyblue")[ind2,]
# }
# 
# if (isFALSE(Silent)) 
#   message("DensityScatter.DDCAL: using ggplot2..")
# DF = data.frame(x = X[ind], y = Y[ind], Colors = cols, 
#                 Cls = cls)
# ggobj = ggplot2::ggplot(DF, ggplot2::aes_string(x = "x", y = "y", 
#                                          col = "Colors")) + ggplot2::geom_point(data = DFnull, 
#                                                                                 mapping = ggplot2::aes_string(x = "x2", y = "y2", col = "Colors"), 
#                                                                                 size = Size, shape = 3, alpha = 0.4) + 
#   ggplot2::geom_point(size = Size, shape = pch, 
#                       alpha = 0.4) + ggplot2::theme(legend.position = "none") + 
#   ggplot2::scale_color_identity()

#taken from plotly variant
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