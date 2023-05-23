DDCAL=function(data, nClusters, minBoundary = 0.1, maxBoundary = 0.45,
               numSimulations = 20, csTolerance = 0.45, csToleranceIncrease = 0.5) {

  # DDCAL
  # DDCAL-Clustering-Algorithm (Density Distribution Cluster Algorithm)
  #
  #  INPUT
  #  data                   [1:n] Numeric vector of the input data
  #  nClusters              Scalar, number of clusters
  #  OPTIONAL
  #  minBoundary            Scalar, in the range (0,1), gives the lower boundary (in percent), for the simulation
  #  maxBoundary            Scalar, in the range (0,1), gives the upper boundary (in percent), for the simulation
  #  csTolerance            Scalar, in the range (0,1). Gives cluster size tolerance factor.
  #                         The necessary cluster size is defined by
  #                         (dataSize/nClusters - dataSize/nClusters * csTolerance)
  #  csToleranceIncrease    Scalar, in the range (0,1), gives the procentual increase of the
  #                         csTolerance-factor, if some clusters did not reach the necessary size.
  #  OUTPUT
  #  labels                 [1:n] Numeric vector, containing the labels for the input data points
  #
  #  Author LB 03/2023
  ##############




  if(!is.vector(data)) {
    if(is.matrix(data)) {
      data = data[,1]   # If a matrix is given as Input, use the first row
      warning("DDCAL: Data was given as a matrix, only first column will be used")
    } else {
      stop("DDCAL: Data has to be either a vector or a matrix")
    }
  }

  dataSizeComplete = length(data)
  ind_na = rep(TRUE, length(data))
  if(isTRUE(any(!is.finite(data)))) {
    warning("DDCAL: Non finite values are found in input data. Non finite values will not be clustered")
    ind_na = is.finite(data)
    data = data[ind_na]
  }



  if(nClusters <= 0) {
    stop("Number of Clusters (nClusters) has to be greater than 0")
  }
  if(numSimulations <= 0) {
    warning("Number of Simulations (numSimulations) has to be greater than 0. Setting numSimulations = 20")
    numSimulations = 20
  }
  if(minBoundary < 0 || minBoundary > 1) {
    warning("minBoundary has to be between 0 and 1. Setting minBoundary = 0.1")
    minBoundary = 0.1
  }
  if(maxBoundary < 0 || maxBoundary > 1) {
    warning("maxBoundary has to be between 0 and 1. Setting maxBoundary = 0.45")
    maxBoundary = 0.45
  }
  if(csTolerance < 0 || csTolerance > 1) {
    warning("csTolerance has to be between 0 and 1. Setting csTolerance = 0.45")
    csTolerance = 0.45
  }
  if(csToleranceIncrease < 0 || csToleranceIncrease > 1) {
    warning("csToleranceIncrease has to be between 0 and 1. Setting csToleranceIncrease = 0.5")
    csToleranceIncrease = 0.5
  }

  # cs = cluster size
  csToleranceTmp = csTolerance
  clustersArray = seq(from = 0, to = nClusters - 1, by = 1)
  dataSize = length(data)
  dataMinIndex = 1
  dataMaxIndex = dataSize
  labels = rep(-1, dataSize)
  average_cs = dataSize/nClusters
  minData = min(data)
  maxData = max(data)

  if(nClusters == 1) {
    return(rep(1, dataSize))
  }

  while(length(clustersArray) > 0 && dataMinIndex < dataMaxIndex) {
    dataTmp  = data[dataMinIndex:dataMaxIndex]
    minDataTmp = min(dataTmp)
    maxDataTmp = max(dataTmp)
    if(minDataTmp != maxDataTmp) {
      normFrequencies = (dataTmp - min(dataTmp)) / (max(dataTmp) - min(dataTmp))
    } else {
      normFrequencies = rep(0, length(dataTmp))
    }
    min_cs = average_cs - average_cs * csToleranceTmp

    foundBestIndexLow = FALSE
    foundBestIndexUp = FALSE

    if(minBoundary == maxBoundary) {
      numSimulations = 1
    }
    boundaries = seq(from=minBoundary, to=maxBoundary, length.out=numSimulations) # Why calculate not outside loop?

    for(b in 1:length(boundaries)) {
      lowerBound = 0
      upperBound = 0
      for(i in 1:length(normFrequencies)) {
        if(normFrequencies[i] <= (boundaries[b])) {
          lowerBound = lowerBound + 1
        }
        if(normFrequencies[i] >= (1 - boundaries[b])) {
          upperBound = upperBound + 1
        }
      }
      # Find potential cluster boundary for the left Cluster
      while(dataMinIndex + lowerBound <= dataMaxIndex && lowerBound > 1) { # if lowerBound is 1, than sd is not defined
        stdDevCurrent = sd(data[dataMinIndex:(dataMinIndex + lowerBound - 1)])
        indexNextGap = getNextGap(data, dataMinIndex + lowerBound, dataMinIndex, dataMaxIndex - 1, 1)
        stdDev_Next = sd(data[dataMinIndex:indexNextGap])
        if(stdDev_Next < stdDevCurrent) {
          lowerBound = lowerBound + indexNextGap - (dataMinIndex + lowerBound)
        } else {
          break # potential cluster boundary for left cluster found
        }
      }
      # Find potential cluster boundary for the right Cluster
      while(dataMaxIndex - upperBound >= dataMinIndex && upperBound > 1) { # if upperBound is 1, than sd is not defined
        stdDevCurrent = sd(data[(dataMaxIndex - upperBound + 1):dataMaxIndex])
        indexNextGap = getNextGap(data, dataMaxIndex - upperBound, dataMinIndex, dataMaxIndex - 1, -1)
        stdDev_Next = sd(data[indexNextGap:dataMaxIndex])
        if(stdDev_Next < stdDevCurrent) {
          upperBound = upperBound + dataMaxIndex - (indexNextGap + upperBound)
        } else {
          break # potential cluster boundary for right cluster found
        }
      }

      if(min_cs > lowerBound || min_cs > upperBound) {
        next # Either left or right cluster is not big enough
      } else {
        if(abs(lowerBound - average_cs) <= abs(upperBound - average_cs)) { # lowerBound <= upperBound should also be working
          foundBestIndexLow = TRUE
        } else {
          foundBestIndexUp = TRUE
        }
        break
      }
    }
    if(!foundBestIndexLow & !foundBestIndexUp) {
      csToleranceTmp = csToleranceTmp + csToleranceTmp * csToleranceIncrease
      next # No cluster of minimum required size was found => give cluster size more tolerance
    }
    if(foundBestIndexLow) {
      labels[dataMinIndex:(dataMinIndex + lowerBound)] = clustersArray[1]
      dataMinIndex = dataMinIndex + lowerBound
      clustersArray = clustersArray[-1]
    } else {
      labels[(dataMaxIndex - upperBound):dataMaxIndex] = clustersArray[length(clustersArray)]
      dataMaxIndex = dataMaxIndex - upperBound
      clustersArray = clustersArray[-length(clustersArray)]
    }
    if(length(clustersArray) == 1) {
      labels[dataMinIndex:dataMaxIndex] = clustersArray[1]
      break
    }
    csToleranceTmp = csTolerance
    average_cs = (dataMaxIndex - dataMinIndex) / length(clustersArray)
  }

  # Labels are for sorted data => unsort the labels before returning
  indicesOfSortedValues = order(data)
  labelsUnsorted = numeric(length(indicesOfSortedValues))
  l = 1
  for(j in 1:length(indicesOfSortedValues)){
    labelsUnsorted[indicesOfSortedValues[j]] = labels[l]
    l = l + 1
  }
  labelsUnsorted_na = rep(NaN, length(ind_na))
  labelsUnsorted_na[ind_na] = labelsUnsorted
  return(labelsUnsorted_na)
}