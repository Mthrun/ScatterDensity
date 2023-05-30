DensityScatter.DDCAL = function (X, Y, xlab, ylab, SDHorPDE = TRUE, Plotter = "native", 
            Silent = FALSE, Marginals = FALSE, pch = 10, Size = 1, BW = TRUE, 
            PDEsample = 5000, lwd = 2, na.rm = TRUE,Polygon, ...) 
  {
    
    if (length(X) != length(Y)) {
      stop("DensityScatter.DDCAL: length of X points does not equal lengt of Y-points")
    }
    if (missing(xlab)) 
      xlab = deparse1(substitute(X))
    if (missing(ylab)) 
      ylab = deparse1(substitute(Y))
    
    if(isTRUE(na.rm)) {
      ind.na = is.finite(X) & is.finite(Y)
      X = X[ind.na]
      Y = Y[ind.na]
    }
    Bool = FALSE
    if (length(X) > PDEsample && isFALSE(SDHorPDE)) {
      Bool = TRUE
      ind = sample(1:length(X), PDEsample)
    }
    else {
      ind = 1:length(X)
    }
    if (isFALSE(Silent)) 
      message("DensityScatter.DDCAL: Estimating Density...")
    if (isTRUE(SDHorPDE)) {
      Dens = SmoothedDensitiesXY(X = X[ind], Y = Y[ind], PlotIt = FALSE)$Densities
    }
    else {
      Dens = PDEscatter(x = X[ind], y = Y[ind], PlotIt = -1)$Densities
    }
    
    if (isFALSE(Silent)) 
      message("DensityScatter.DDCAL: Estimating colors...")
    cls = DDCAL(Dens, nClusters = 12, minBoundary = 0.2, 
                maxBoundary = 0.6, numSimulations = 20, 
                csTolerance = 0.45, csToleranceIncrease = 0.5)
    xlim = c(min(X, na.rm = T), max(X, na.rm = T))
    ylim = c(min(Y, na.rm = T), max(Y, na.rm = T))
    cols <- rep(NA_character_, length(ind))
    ncolors = 1
    ncolors2 = 1
    colpalette = colorRampPalette(c("navyblue", "darkblue", rep("blue", 
                                  ncolors2), rep("turquoise", ncolors2), rep("green", ncolors), 
                                  rep("chartreuse", ncolors), rep("yellow", ncolors), rep("orange", 
                                  ncolors), rep("red", ncolors), rep("darkred", ncolors2)))
    
    colpal <- cut(cls, length(cls), labels = FALSE)
    cols <- colpalette(length(cls))[colpal]
    if (Plotter == "native") {
      if (isTRUE(Marginals)) {
        def.par <- par(no.readonly = TRUE)
        on.exit(par(def.par))
        
        m <- graphics::layout(matrix(c(1, 2, 2, 2, 1, 2, 
                                       2, 2, 1, 2, 2, 2, 4, 3, 3, 3), 4, 4))
        par(fig = c(0, 0.8, 0, 0.8), 
            mar = c(5, 4, 0, 0))
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
             cex = Size, ...)
        points(x = X[ind], y = Y[ind], col = cols, pch = pch, 
               cex = Size, ...)
      }
      else {
        plot(x = X, y = Y, xlab = xlab, ylab = ylab, xlim = xlim, 
             ylim = ylim, col = cols, pch = pch, cex = Size, ...)
      }
      if(!missing(Polygon)){
        points(Polygon[,1],Polygon[,2],lwd=3,col="magenta",type="l")
      }
      
      if (isTRUE(Marginals)) {
        par(fig = c(0, 0.8, 0.8, 
                    1), mar = c(0, 4, 4, 0), new = TRUE)
        V = DataVisualizations::ParetoDensityEstimation(X[ind], PlotIt = F)
        plot(V$kernels, V$paretoDensity, xlim = xlim, ylim = c(0, max(V$paretoDensity)), 
             ylab = "PDE", type = "l", lwd = lwd, axes = FALSE, bty = "n")
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
    }
    else if (Plotter == "ggplot2") {
      if (isFALSE(Silent)) 
        message("DensityScatter.DDCAL: Preparing for ggplot2..")
      if (isTRUE(Bool)) {
        if (length(X) > 2 * PDEsample) {
          mod = ClusterR::MiniBatchKmeans(cbind(X, Y), 
                                          2 * PDEsample)
          Centroids = mod$centroids
          DFnull = data.frame(x2 = Centroids[, 1], y2 = Centroids[, 
                                                                  2], Colors = "navyblue")
        }
        else {
          ind2 = 1:length(X)
          DFnull = data.frame(x2 = X, y2 = Y, Colors = "navyblue")[ind2, 
          ]
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
                                                 group = "Cls", col = "Colors"), alpha = 0.05) + 
          ggplot2::geom_point(size = Size, shape = pch) + 
          ggplot2::theme(legend.position = "none") + ggplot2::scale_color_identity()
      }
      if (isTRUE(BW)) 
        ggobj = ggobj + ggplot2::theme_bw()
      if (isTRUE(Marginals)) {
        ggobj2 <- ggExtra::ggMarginal(ggobj, type = "density")
        #print(ggobj2)
        return(ggobj2)
      }
      else {
        #print(ggobj)
        return(ggobj)
      }
    }
    else {
      if (isFALSE(Silent)) 
        message("DensityScatter.DDCAL: plotly not implemented")
    }
  }
