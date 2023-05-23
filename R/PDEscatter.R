PDEscatter <-
function (x, y, SampleSize, na.rm = FALSE, PlotIt = TRUE, ParetoRadius, 
    sampleParetoRadius, NrOfContourLines = 20, Plotter = "native", 
    DrawTopView = TRUE, xlab = "X", ylab = "Y", main = "PDEscatter", 
    xlim, ylim, Legendlab_ggplot = "value") 
{
    x = checkFeature(x, varname = "x", Funname = "PDEscatter")
    y = checkFeature(y, varname = "y", Funname = "PDEscatter")
    if (identical(x, y)) {
        stop("PDEscatter: Variable x is identical to variable y. Please check input.")
    }
    isnumber = function(x) return(is.numeric(x) & length(x) == 
        1)
    if (missing(ParetoRadius)) {
        ParetoRadius = 0
    }
    if (is.null(ParetoRadius)) {
        ParetoRadius = 0
    }
    if (!isnumber(ParetoRadius)) 
        stop("PDEscatter: \"ParetoRadius\" is not a numeric number of length 1. Please change Input.")
    if (missing(SampleSize)) {
        SampleSize = -1
    }
    if (missing(sampleParetoRadius)) {
        sampleParetoRadius = round(sqrt(500000000), -3)
    }
    if (!isnumber(SampleSize)) 
        stop("PDEscatter: \"SampleSize\" is not a numeric number of length 1. Please change Input.")
    if (!isnumber(sampleParetoRadius)) 
        stop("PDEscatter: \"sampleParetoRadius\" is not a numeric number of length 1. Please change Input.")
    if (!isnumber(NrOfContourLines)) 
        stop("PDEscatter: \"NrOfContourLines\" is not a numeric number of length 1. Please change Input.")
    if (!is.logical(na.rm)) 
        stop("PDEscatter: \"na.rm\" is expected to be either TRUE or FALSE")
    if (!is.logical(PlotIt)) {
        if (!(PlotIt == -1)) 
            stop("PDEscatter: \"PlotIt\" is expected to be either TRUE, FALSE or -1.")
    }
    if (!is.logical(DrawTopView)) 
        stop("PDEscatter: \"DrawTopView\" is expected to be either TRUE or FALSE")
    toRange = function(data, lower, upper) {
        data <- as.matrix(data)
        if (lower == upper) {
            stop("interval width can not be 0!")
        }
        if (lower > upper) {
            temp <- upper
            upper <- lower
            lower <- upper
        }
        range <- upper - lower
        n <- dim(data)[1]
        d <- dim(data)[2]
        if ((n == 1) & (d > 1)) {
            data <- t(data)
            wasRowVector <- 1
        }
        else {
            wasRowVector <- 0
        }
        nRow <- dim(data)[1]
        nCol <- dim(data)[2]
        min <- apply(data, 2, min, na.rm = TRUE)
        min <- matrix(min, nRow, nCol, byrow = TRUE)
        max <- apply(data, 2, max, na.rm = TRUE)
        max <- matrix(max, nRow, nCol, byrow = TRUE)
        range <- max - min
        range[range == 0] <- 1
        scaleData <- (data - min)/range
        scaleData <- lower + scaleData * (upper - lower)
        if (wasRowVector == 1) {
            scaleData = t(scaleData)
        }
        return(scaleData)
    }
    if (isTRUE(na.rm)) {
        noNaNInd <- which(is.finite(x) & is.finite(y))
        x = x[noNaNInd]
        y = y[noNaNInd]
    }
    nData <- length(x)
    if (SampleSize > 0) {
        if (SampleSize < nData) {
            sampleInd = sample(1:nData, size = SampleSize)
            x = x[sampleInd]
            y = y[sampleInd]
        }
    }
    if (missing(xlim)) 
        xlim = c(min(x, na.rm = T), max(x, na.rm = T))
    if (missing(ylim)) 
        ylim = c(min(y, na.rm = T), max(y, na.rm = T))
    data <- cbind(x, y)
    percentdata <- toRange(data, 0, 100)
    nData <- length(x)
    if (sampleParetoRadius < nData) {
        par_sampleInd = sample(1:nData, size = sampleParetoRadius, 
            replace = FALSE)
        sampleData4radius = percentdata[par_sampleInd, ]
    }
    else {
        sampleData4radius = percentdata
    }
    if (!requireNamespace("parallelDist", quietly = TRUE)) {
        message("Subordinate package (parallelDist) is missing. No computations are performed.\nPlease install the package which is defined in \"Suggests\". Falling back to dist().")
        DataDists = as.matrix(dist(sampleData4radius, method = "euclidean", 
            diag = TRUE))
        Dists = as.vector(DataDists[upper.tri(DataDists, diag = F)])
    }
    else {
        Dists = parallelDist::parDist(sampleData4radius, method = "euclidean", 
            diag = F, upper = F)
        Dists = as.vector(Dists)
    }
    if (ParetoRadius <= 0) {
        if (nData < 500) {
            ParetoRadius <- quantile(Dists, 6/100, type = 5, 
                na.rm = TRUE)
        }
        else {
            ParetoRadius <- c_quantile(Dists[is.finite(Dists)], 
                6/100)
        }
        if (ParetoRadius == 0) {
            if (nData < 500) {
                ParetoRadius <- quantile(Dists, 20/100, type = 5, 
                  na.rm = TRUE)
            }
            else {
                ParetoRadius <- c_quantile(Dists[is.finite(Dists)], 
                  20/100)
            }
            if (ParetoRadius == 0) {
                stop(paste0("Estimation of Radius(", ParetoRadius, 
                  ") for two-dimensional density not possible. Please provide ParetoRadius manually."))
            }
            else {
                warning(paste0("Estimation of Radius(", ParetoRadius, 
                  ") for two-dimensional density may not work properly. You can provide ParetoRadius manually."))
            }
        }
    }
    inPSpheres = inPSphere2D(percentdata, ParetoRadius)
    Matrix3D = cbind(x, y, inPSpheres)
    if (PlotIt == -1) 
        return(list(X = x, Y = y, Densities = inPSpheres, Matrix3D = Matrix3D, 
            ParetoRadius = ParetoRadius, Handle = NULL))
    if(requireNamespace("DataVisualizations")){
      plt = DataVisualizations::zplot(x = x, y = y, z = inPSpheres, DrawTopView, NrOfContourLines, 
                                      TwoDplotter = Plotter, xlim = xlim, ylim = ylim)
    }else{
      return(list(X = x, Y = y, Densities = inPSpheres, Matrix3D = Matrix3D, 
                  ParetoRadius = ParetoRadius, Handle = NULL))
    }

    if (mode(plt) == "character") {
        return(list(X = x, Y = y, Densities = inPSpheres, Matrix3D = Matrix3D, 
            ParetoRadius = ParetoRadius, Handle = NULL))
    }
    if (DrawTopView) {
        switch(Plotter, ggplot = {
          if(requireNamespace("ggplot2")){
            plt <- plt + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::labs(title = main, 
                fill = Legendlab_ggplot) + ggplot2::theme(panel.background = ggplot2::element_blank())
            if (isTRUE(PlotIt)) print(plt)
          }else{
            plt=NULL
          }
        }, native = {
            title(main = main, xlab = xlab, ylab = ylab)
            plt <- "Native does not have a Handle"
            if (!isTRUE(PlotIt)) warning("for native plotting cannot be disabled")
        }, plotly = {
            if(requireNamespace("plotly")){
            plt <- plotly::layout(plt, xaxis = list(title = xlab), 
                yaxis = list(title = ylab), title = main)
            }else{
              plt=NULL
            }
            if (isTRUE(PlotIt)) print(plt)
        })
    }
    else {
        switch(Plotter, ggplot = {
            print("Plotly plot is used because ggplot is not implemented for option DrawTopView=FALSE.")
            if(requireNamespace("plotly")){
            plt <- plotly::layout(plt, scene = list(xaxis = list(title = xlab), 
                yaxis = list(title = ylab), zaxis = list(title = "PDE"), 
                title = main))
            if (isTRUE(PlotIt)) print(plt)
            }else{
              plt=NULL
            }
     
        }, native = {
            print("Plotly plot is used because native is not implemented for option DrawTopView=FALSE.")
            requireNamespace("plotly")
            plt <- plotly::layout(plt, scene = list(xaxis = list(title = xlab), 
                yaxis = list(title = ylab), zaxis = list(title = "PDE"), 
                title = main))
            if (isTRUE(PlotIt)) print(plt)
        }, plotly = {
            if(requireNamespace("plotly")){
            plt <- plotly::layout(plt, scene = list(xaxis = list(title = xlab), 
                yaxis = list(title = ylab), zaxis = list(title = "PDE"), 
                title = main))
            if (isTRUE(PlotIt)) print(plt)
            }else{
              plt=NULL
            }
        })
    }
    return(invisible(list(X = x, Y = y, Densities = inPSpheres, 
        Matrix3D = Matrix3D, ParetoRadius = ParetoRadius, Handle = plt)))
}
