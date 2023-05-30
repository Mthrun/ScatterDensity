inPSphere2D = function (data, paretoRadius = NULL) 
{
    if (is.null(paretoRadius)){
      if(requireNamespace("DataVisualizations")){
        paretoRadius <- DataVisualizations::ParetoRadius(data)
      } else{
        warning("inPSphere2D: Please install package 'DataVisualizations' or set paretoRadius.")
        paretoRadius=1
      }
    }
      
    nData <- nrow(data)
    nVar <- ncol(data)
    x <- data[, 1]
    y <- data[, 2]
    noNaNInd <- (!is.nan(x) & !is.nan(y))
    x <- x[noNaNInd]
    y <- y[noNaNInd]
    xMin <- min(x)
    xMax <- max(x)
    yMin <- min(y)
    yMax <- max(y)
    xBinWidth <- paretoRadius
    yBinWidth <- paretoRadius
    xedge <- seq(xMin, (xMax + xBinWidth), by = xBinWidth)
    yedge <- seq(yMin, (yMax + yBinWidth), by = yBinWidth)
    if (min(x) < xedge[1]) 
        xedge <- rbind(min(x), xedge)
    if (max(x) > xedge[length(xedge)]) 
        xedge <- rbind(xedge, max(x))
    if (min(y) < yedge[1]) 
        yedge <- rbind(min(y), yedge)
    if (max(y) > yedge[length(yedge)]) 
        yedge <- rbind(yedge, max(y))
    e <- hist(x, xedge, plot = FALSE)
    xBinNr <- findInterval(x, e$breaks)
    e <- hist(y, yedge, plot = FALSE)
    yBinNr <- findInterval(y, e$breaks)
    nrXBins <- length(xedge)
    nrYBins <- length(yedge)
    nInPsphere = c_inPSphere2D(data, xBinNr, yBinNr, nrXBins, 
        nrYBins, nData, paretoRadius)
    return(nInPsphere)
}
