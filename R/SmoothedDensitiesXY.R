SmoothedDensitiesXY = function (X, Y, nbins, lambda, Xkernels, Ykernels, PlotIt = FALSE) 
{
    if (!requireNamespace("pracma")) 
        stop("pracma package is missing")
    if (missing(nbins)) {
        nbins = c(min(length(unique(X)), 2000), min(length(unique(Y)), 
            2000))
    }
    else if (is.null(nbins)) {
        nbins = c(min(length(unique(X)), 2000), min(length(unique(Y)), 
            2000))
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
    }
    else if (missing(Xkernels)) {
        stop("SmoothedDensitiesXY: Ykernels given but Xkernels missing, case not implemented!")
    }
    else if (missing(Ykernels)) {
        stop("SmoothedDensitiesXY: Xkernels given but Ykernels missing, case not implemented!")
    }
    else {
        minx = min(X)
        maxx = max(X)
        miny = min(Y)
        maxy = max(Y)
        edges1 = seq(from = minx, to = maxx, length.out = nbins[1] + 
            1)
        edges1 = c(-Inf, edges1[2:(length(edges1) - 1)], Inf)
        edges2 = seq(from = miny, to = maxy, length.out = nbins[2] + 
            1)
        edges2 = c(-Inf, edges2[2:(length(edges2) - 1)], Inf)
    }
    bin = matrix(0, n, 2)
    bin[, 2] = pracma::histc(X, edges1)$bin
    bin[, 1] = pracma::histc(Y, edges2)$bin
    H = pracma::accumarray(bin, rep(1, nrow(bin)), nbins[c(2, 
        1)])/n
    G = smooth1D(H, nbins[2]/lambda)
    hist_F_2D = t(smooth1D(t(G), nbins[1]/lambda))
    MaxF = max(as.vector(hist_F_2D))
    hist_F_2D = hist_F_2D/MaxF
    m = dim(hist_F_2D)[1]
    ind = (bin[, 2] - 1) * m + bin[, 1]
    XYDensity = hist_F_2D[ind]
    Density = matrix(0, OrigN, 1)
    Density[NoNaNInd] = XYDensity
    Density = as.vector(Density)
    if (isTRUE(PlotIt)&requireNamespace("DataVisualizations")) {
        DataVisualizations::zplot(X, Y, Density)
    }
    return(list(Densities = as.vector(Density), Xkernels = Xkernels, 
        Ykernels = Ykernels, hist_F_2D = hist_F_2D, ind = ind))
}
