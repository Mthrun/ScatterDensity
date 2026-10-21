SmoothedDensitiesXY = function (X, Y, nbins, lambda, Xkernels, Ykernels, isXDiscrete=FALSE, isYDiscrete=FALSE, Compute="Cpp", PlotIt = FALSE) 
{
    if (!requireNamespace("pracma")) 
        stop("pracma package is missing")
    if (missing(nbins)) {
        nbins = c(min(length(unique(X)), 250), min(length(unique(Y)), 
                 250))
    }
    else if (is.null(nbins)) {
        nbins = c(min(length(unique(X)), 250), min(length(unique(Y)), 
                  250))
    }
    else if (length(nbins) == 1) {
        nbins = c(nbins, nbins)
    }
    else {
    }
    if (missing(lambda)) {
        lambda = 20
    }
    else if (is.null(lambda)) {
        lambda = 20
    }
    else {
    }
  
  Compute=tolower(Compute)
  Compute=gsub("\\+","p",Compute)
  
    X = as.vector(X)
    Y = as.vector(Y)
    OrigN = length(X)
    NoNaNInd <- which(is.finite(X) & is.finite(Y))
    X = X[NoNaNInd]
    Y = Y[NoNaNInd]
    n = length(X)
    ny = length(Y)
    if (n != ny) 
        stop(" SmoothedDensitiesXY: length(X) is not the same as length(Y)")
    
    if (!isXDiscrete && !isYDiscrete) {
      if (missing(Xkernels) & missing(Ykernels)) {
        minx = min(X)
        maxx = max(X)
        miny = min(Y)
        maxy = max(Y)
        edges1 = seq(from = minx, to = maxx, length.out = nbins[1] + 
                       1)
        Xkernels = edges1[1:(length(edges1) - 1)] + 0.5 * diff(edges1)
        edges1 = c(-Inf, edges1[2:(length(edges1) - 1)], Inf)
        edges2 = seq(from = miny, to = maxy, length.out = nbins[2] + 
                       1)
        Ykernels = edges2[1:(length(edges2) - 1)] + 0.5 * diff(edges2)
        edges2 = c(-Inf, edges2[2:(length(edges2) - 1)], Inf)
        
        bin = matrix(0, n, 2)
        
        bin[, 2] = pracma::histc(X, edges1)$bin
        bin[, 1] = pracma::histc(Y, edges2)$bin
        
        switch (Compute,
                cpp = {
                  H=accumarray_rcpp(bin, rep(1, nrow(bin)), nbins[c(2, 
                                                                    1)])/n
                  G = smooth1D_C(H, nbins[2]/lambda)
                  hist_F_2D = t(smooth1D_C(t(G), nbins[1]/lambda))
                },
                r = {
                  H = pracma::accumarray(bin, rep(1, nrow(bin)), nbins[c(2, 
                                                                         1)])/n
                  G = smooth1D(H, nbins[2]/lambda)
                  hist_F_2D = t(smooth1D(t(G), nbins[1]/lambda))
                },
                parallel = {
                  H=accumarray_rcpp(bin, rep(1, nrow(bin)), nbins[c(2, 
                                                                    1)])/n
                  G = smooth1D_parallel(H, nbins[2]/lambda, nbins[2],smooth_rows=TRUE)
                  hist_F_2D = smooth1D_parallel(G,nbins[1]/lambda, nbins[1],smooth_rows=FALSE)
                  #hist_F_2D = t(smooth1D_parallel(t(G), nbins[1]/lambda))
                },
                {
                  stop("SmoothedDensitiesXY:Incorrect Compute parameter selected. Options are 'R','Cpp,'Parallel'")
                }
        )
        
        #MaxF = max(as.vector(hist_F_2D))
        #hist_F_2D = hist_F_2D/MaxF
        dx = diff(Xkernels)
        dy = diff(Ykernels)
        int = sum(hist_F_2D[1:nrow(hist_F_2D)-1, 1:ncol(hist_F_2D)-1] * outer(dx, dy))
        hist_F_2D = hist_F_2D / int
        
        m = dim(hist_F_2D)[1]
        ind = (bin[, 2] - 1) * m + bin[, 1]
        XYDensity = hist_F_2D[ind]
        Density = matrix(0, OrigN, 1)
        Density[NoNaNInd] = XYDensity
        Density = as.vector(Density)
        
      }#both kernels are not missing
      else if (missing(Xkernels)) {
        stop("SmoothedDensitiesXY: Ykernels given but Xkernels missing, case not implemented!")
      }
      else if (missing(Ykernels)) {
        stop("SmoothedDensitiesXY: Xkernels given but Ykernels missing, case not implemented!")
      }
      else {#both kernels are given
        nbinsX = length(Xkernels)
        nbinsY = length(Ykernels)
        
        nbins <- c(nbinsX, nbinsY)
        #kernels sind die zentrien
        #mit diff komme ich zu den grenzen
        # -inf und inf stell sicher dass ich alle datenpunkte in einem bin habe
        edges1 <- c(-Inf, Xkernels[1:(nbinsX-1)]+0.5*diff(Xkernels), Inf)
        edges2 <- c(-Inf, Ykernels[1:(nbinsY-1)]+0.5*diff(Ykernels), Inf)
        #3 kernels gibt 2 edges sowie die beiden ausengrenzen, also 3 bins
        #4 kernels gibt 3 edges sowie die beiden ausengrenzen, also 4 bins
        
        # binning
        bin <- matrix(0, nrow = n, ncol = 2)
        bin[, 2] <- findInterval(X, vec = edges1, all.inside = TRUE)
        bin[, 1] <- findInterval(Y, vec = edges2, all.inside = TRUE)
        
        # Histogramm bauen
        H <- matrix(0, nrow = nbins[2], ncol = nbins[1])
        for (i in 1:n) {
          if (!is.na(bin[i, 1]) && !is.na(bin[i, 2]) && bin[i, 1] > 0 && bin[i, 2] > 0) {
            H[bin[i, 2], bin[i, 1]] <- H[bin[i, 2], bin[i, 1]] + 1
          }
        }
        H <- H / n
        
        # Glättung
        switch (Compute,
                cpp = {
                  G = smooth1D_C(H, nbins[2]/lambda)
                  hist_F_2D = t(smooth1D_C(t(G), nbins[1]/lambda))
                },
                r = {
                  
                  G <- smooth1D(H, nbins[2] / lambda)
                  hist_F_2D <- t(smooth1D(t(G), nbins[1] / lambda))
                },
                parallel = {
                  G = smooth1D_parallel(H, nbins[2]/lambda, nbins[2],smooth_rows=TRUE)
                  hist_F_2D = smooth1D_parallel(G,nbins[1]/lambda, nbins[1],smooth_rows=FALSE)
                },
                {
                  stop("SmoothedDensitiesXY:Incorrect Compute parameter selected. Options are 'R','Cpp,'Parallel'")
                }
        )  
        
        # Normierung:
        #hist_F_2D <- hist_F_2D / max(hist_F_2D, na.rm = TRUE)
        dx = diff(Xkernels)
        dy = diff(Ykernels)
        int = sum(hist_F_2D[1:nrow(hist_F_2D)-1, 1:ncol(hist_F_2D)-1] * outer(dx, dy))
        hist_F_2D = hist_F_2D / int
        
        # Dichte an jedem Punkt
        ind <- (bin[, 1] - 1) * nbins[2] + bin[, 2]
        Density <- hist_F_2D[ind]
        
      }# end if (missing(Xkernels) & missing(Ykernels)) 
      
      
      if (isTRUE(PlotIt)&requireNamespace("DataVisualizations")) {
        DataVisualizations::zplot(X, Y, Density)
      }
      return(list(Densities = as.vector(Density), Xkernels = Xkernels, 
                  Ykernels = Ykernels, GridDensity = hist_F_2D, Points2GridInd = ind))
    } else {
      if (isXDiscrete && isYDiscrete) {
        # X is discrete, Y is discrete
        edges_x=sort(unique(X),decreasing = F)
        edges_y=sort(unique(Y),decreasing = F)
        joint_table=fast_table_num(X,Y,edges_x,edges_y,redefine=F,na.rm = TRUE)
        Density = joint_table / sum(joint_table)
        return(list(
          Densities = as.vector(Density),
          Xkernels = edges_x,
          Ykernels = edges_y,
          GridDensity = Density,
          Points2GridInd = seq_along(Density)
        ))
      } else if (isXDiscrete) {
        # X is discrete, Y is continuous
        t = X
        X = Y
        Y = t
        
        Xkernels = NULL
        if (!missing(Ykernels)) {
          Xkernels = Ykernels
        }
      }
      
      # X continuous, Y discrete
      if (missing(Xkernels) || is.null(Xkernels)) {
        minx = min(X)
        maxx = max(X)
        edges1 = seq(from = minx, to = maxx, length.out = nbins[1] + 
                       1)
        Xkernels = edges1[1:(length(edges1) - 1)] + 0.5 * diff(edges1)
        edges1 = c(-Inf, edges1[2:(length(edges1) - 1)], Inf)
      } else {
        nbinsX = length(Xkernels)
        
        edges1 <- c(-Inf, Xkernels[1:(nbinsX-1)]+0.5*diff(Xkernels), Inf)
      }
      edges2=sort(unique(Y),decreasing = F)
      edges2=c(-Inf, (head(edges2, -1) + tail(edges2, -1)) / 2, Inf)
      
      #joint_table=Tiles::fast_table_num(X,Y,edges1,redefine=F,na.rm = TRUE)
      # since currently the output is transposed:
      joint_table=t(fast_table_num(X,Y,edges1,edges2,redefine=F,na.rm = TRUE))
      
      joint_table <- joint_table / n
      switch (Compute,
              cpp = {
                hist_F_2D <- smooth1D_C(joint_table, nbins[1] / lambda)
              },
              r = {
                hist_F_2D <- smooth1D(joint_table, nbins[1] / lambda)
              },
              parallel = {
                stop("SmoothedDensitiesXY: Compute Parameter 'Parallel' is not supported for discrete variables")
                #hist_F_2D <- smooth1D_parallel(joint_table, nbins[1] / lambda, length(edges2), smooth_rows=TRUE)
              },
              {
                stop("SmoothedDensitiesXY:Incorrect Compute parameter selected. Options are 'R','Cpp,'Parallel'")
              }
      )  
      
      # Normierung:
      #hist_F_2D <- hist_F_2D / max(hist_F_2D, na.rm = TRUE)
      dx = diff(Xkernels)
      int = sum(colSums(sweep(hist_F_2D[1:nrow(hist_F_2D)-1,], 1, dx, `*`)))
      hist_F_2D = hist_F_2D / int
      
      
      
      if (isXDiscrete) {
        # reflip the output
        hist_F_2D = t(hist_F_2D)
        
        Ykernels = Xkernels
        Xkernels = edges2
      } else {
        Ykernels = edges2
      }
      return(list(
        Densities = as.vector(hist_F_2D),
        Xkernels = Xkernels,
        Ykernels = Ykernels,
        GridDensity = hist_F_2D,
        Points2GridInd = seq_along(hist_F_2D)
      ))
    }
    
    
}
