UniquePoints <- function(Datapoints, Eps) {
  # 
  # V <- UniquePoints(Datapoints)
  # return only the unique points in Datapoints
  # 
  # DESCRIPTION
  # This function will return the set of unique points within an epsilon 
  # neighborhood.
  # 
  # INPUT
  # Datapoints[1:n,1:d]   [1:n,1:d] matrix of Datapoints points of dimension d
  #				                the points are in the  rows
  #
  # Eps                   Optional,scalarabove zero that defines minimum non-
  #                       identical euclidean distance between two points
  # 
  # OUTPUT
  # a list V containing:
  # Unique[1:k,1:d]        the Datapoints points  without duplicate points
  # IsDuplicate[1:n]       is the datapoint a unique or a doublet
  # (IsDuplicate[1:n,1:n]   for i!=j IsDuplicate[i,j]== 1  if Datapoints[i,] == Datapoints[j,]    IsDuplicate[i,i]==0)
  # UniqueInd[1:k]		     an index vector 1:k such that
  #                        Unique == Datapoints[UniqueInd,], it has k
  #                        non-consecutive numbers or labels, each label defines
  #                        a row number within Datapoints[1:n,1:d] of a unique
  #                        data point
  # Uniq2DatapointsInd[1:n]    an index vector 1:n. It has k unique index
  #                            numbers representing the arbitrary labels. Each
  #                            labels is mapped uniquely to a point in Unique.
  #                            Logically in a way such that
  #                            Datapoints == Unique[Uniq2DatapointsInd,]
  #                            (will not work directly in R this way) 
  #
  # Description: Euclidean distance is computed and used within. Setting 
  # \code{Eps} to a very small number results in the identification of unique
  # data points. Setting epsilon to a higher number results in the definition of
  # mesh points within an d-dimensional R-ball grap
  # 
  # DETAILS
  # Careful! This function is mathematically not uniquely defined! There can be
  # different outcomes if data is shuffled. This is due to the order of the
  # data, which determines the rank in distances, where the first occurence is
  # always regarded as the unique point, and the subsequent close neighbors are
  # doublets!
  # 
  #MT 06/2022

  if(missing(Eps)){                 # ab dieser Distanz zwischen 2 punkten sind diese identisch
    Eps = 0.00001
  }
  #if(inherits(Datapoints,"matrix")) {
  if(is.vector(Datapoints)) {
      # If Datapoints is a vector.
    Datapoints <- as.matrix(Datapoints)
  }
  AnzPoints = nrow(Datapoints)               # soviele punkte in den Daten
  rownames(Datapoints) = c(1:AnzPoints)      # gib den Zeilen als namen ihren zeilennummer
  
  # Compare each point only with the unique points retained so far.
  # This prevents epsilon chains through points which are themselves duplicates.
  # Distances are calculated on demand, avoiding a complete n x n distance matrix.
  DataMatrix = as.matrix(Datapoints)
  
  IsDuplicate = rep(FALSE,AnzPoints)
  UniqueInd = integer(AnzPoints)
  Uniq2DatapointsInd = integer(AnzPoints)
  NumberOfUniquePoints = 0L
  
  for(i in seq_len(AnzPoints)){
    if(NumberOfUniquePoints == 0L){
      NumberOfUniquePoints = 1L
      UniqueInd[NumberOfUniquePoints] = i
      Uniq2DatapointsInd[i] = NumberOfUniquePoints
    }else{
      CurrentUniqueInd = UniqueInd[seq_len(NumberOfUniquePoints)]
      Differences = sweep(DataMatrix[CurrentUniqueInd,,drop = FALSE],
                          2,
                          DataMatrix[i,],
                          "-")
      SquaredDistances = rowSums(Differences * Differences)
      NearestUnique = which.min(SquaredDistances)
      
      if(sqrt(SquaredDistances[NearestUnique]) < Eps){
        IsDuplicate[i] = TRUE
        Uniq2DatapointsInd[i] = NearestUnique
      }else{
        NumberOfUniquePoints = NumberOfUniquePoints + 1L
        UniqueInd[NumberOfUniquePoints] = i
        Uniq2DatapointsInd[i] = NumberOfUniquePoints
      }
    }
  }
  
  if(!any(IsDuplicate)){    # keine duplikate gefunden
    return(list("Unique"             = Datapoints,
                "UniqueInd"          = c(1:AnzPoints), #sortind
                "Uniq2DatapointsInd" = c(1:AnzPoints),
                "IsDuplicate"        = IsDuplicate))
  }
  
  UniqueInd = UniqueInd[seq_len(NumberOfUniquePoints)]
  UniqueDatapoints = Datapoints[UniqueInd,]
  
  return(list("Unique"             = UniqueDatapoints,
              "UniqueInd"          = UniqueInd,
              "Uniq2DatapointsInd" = Uniq2DatapointsInd,
              "IsDuplicate"        = IsDuplicate))
}
