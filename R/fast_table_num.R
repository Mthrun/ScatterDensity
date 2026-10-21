fast_table_num=function(x,y,edges_x,edges_y,redefine=TRUE,byrow=FALSE,all.inside = FALSE,rightmost.closed=FALSE,sort=FALSE,na.rm=FALSE,names=FALSE,extendOutput=FALSE){
  #tab=fast_table_num(x,y,redefine=F,sort=T)
  #INPUT
  #x,y                  numerical vectors of same length 1:n
  #OPTIONAL
  #edges_x,edges_y      [1:(k_1+1)]  and [1:(k_2+1)]  numerical  vectors defining the specific borders of x and y, default unique(x), unique(x) for categorical scale
  #redefine             TRUE: resets counts in y direction in order from 1:k_2 to k_2:1
  #byrow                logical. If FALSE (the default) the count matrix is filled by columns, otherwise the matrix is filled by rows.
  #all.inside	          boolean; if true, the returned indices are coerced into 1,...,N-1, i.e., 0 is mapped to 1 and N to N-1
  #rightmost.closed	    boolean; if true, the rightmost interval, vec[N-1] .. vec[N] is treated as closed
  #sort                 boolean, if true, edges_x,edges_y are sorted non-decreasingly 
  #na.rm                boolean, if true, only complete observations are taken into account
  #names                boolean, if true, matrix is named by edges_x[1:k_1] and edges_y[1:k_2] (left-side)
  #extendOutput         boolean, if true gives a list of several outputs otherwise only counts_matrix is returned
  #OUTPUT
  # if extendOutput==FALSE
  #counts_matrix            [1:k_2,1:k_1] numeircal matrix of counts
  #if extendOutput==TRUE
  # list V
  #counts_matrix            [1:k_2,1:k_1] numeircal matrix of counts
  # x_idx, y_idx            numerical vectors [1:k_1] and [1:k_2] of counts for x and y based on edges_x,edges_y 
  #
  #author MCT, 05/2025
  #details: edges_x and edges_y must be sorted non-decreasingly, if not given, set sort=T, edges_x,edges_y can contain Inf,-Inf borders. Beware: kernels are centers of bins, edges_x,edges_y are borders of bins
  #example
  #Cls=Hetap$Cls
  #Cls1= Cls+1
  #fast_table_num(Cls,Cls1,redefine = F)==as.matrix(table(Cls,Cls1))
  
  if (isTRUE(na.rm)) {
    bb = is.finite(x) & is.finite(y)
    x = x[bb]
    y = y[bb]
  }
  
  if (missing(edges_x)) {
    edges_x = unique(x)
    nbinsX = length(edges_x)#kategorieller fall
  } else{
    nbinsX = NULL
  }
  if (missing(edges_y)) {
    edges_y = unique(y)
    nbinsY = length(edges_y)#kategorieller fall
  } else{
    nbinsY = NULL
  }
  
  if (isTRUE(sort)) {
    edges_x = sort(edges_x, decreasing = F, na.last = NA)
    edges_y = sort(edges_y, decreasing = F, na.last = NA)
  }
  
  x_idx = findInterval(
    x,
    vec = edges_x,
    all.inside = all.inside,
    rightmost.closed = rightmost.closed,
    left.open = F
  )
  
  y_idx = findInterval(
    y,
    vec = edges_y,
    all.inside = all.inside,
    rightmost.closed = rightmost.closed,
    left.open = F
  )
  
  #wenn user edges setzt umschließen diese bins
  #bins sind damit immer einer weniger als edges, da edges an den grenzen liegen
  if (is.null(nbinsX))
    #
    nbinsX = length(edges_x) - 1
  
  if (is.null(nbinsY))
    nbinsY = length(edges_y) - 1
  
  if (length(unique(x_idx)) > nbinsX) {
    warning("fast_table_num: x counts outside of bins.")
  }
  if (length(unique(y_idx)) > nbinsY) {
    warning("fast_table_num: y counts outside of bins.")
  }
  #complex und einen tick langsamer
  # if(isFALSE(all.inside)){ #edges
  #   if(isFALSE(rightmost.closed)){
  #     if(tail(edges_x,1)!=Inf){
  #       edges_x2=c(edges_x,Inf)
  #     }else{
  #       edges_x2=edges_x
  #     }
  #     if(tail(edges_y,1)!=Inf){
  #       edges_y2=c(edges_y,Inf)
  #     }else{
  #       edges_y2=edges_y
  #     }
  #   }else{#left-most closed, right-,most is open
  #     if(head(edges_x,1)!= -Inf){
  #       edges_x2=c(-Inf,edges_x)
  #     }else{
  #       edges_x2=edges_x
  #     }
  #     if(head(edges_y,1)!= -Inf){
  #       edges_y2=c(-Inf,edges_y)
  #     }else{
  #       edges_y2=edges_y
  #     }
  #   }
  #   x_idx = .bincode(x, breaks = edges_x2, include.lowest = T, right = rightmost.closed)
  #   y_idx = .bincode(y, breaks = edges_y2, include.lowest = T, right = rightmost.closed)
  #
  #   #rightmost.closed=T:below lowest inteveral
  #   #rightmost.closed=F:above highest inteveral
  #   y_idx[is.na(y_idx)] <- 0
  #   x_idx[is.na(x_idx)] <- 0
  # }else{#kernels
  #   x_idx = .bincode(x, breaks = edges_x, include.lowest = T, right = rightmost.closed)
  #   y_idx = .bincode(y, breaks = edges_y, include.lowest = T, right = rightmost.closed)
  # }
  
  #print(x_idx)
  #print(y_idx)
  #by renaming makes sure that the later ordering by column is correct
  if (isTRUE(redefine))
    y_idx = FCPS::ClusterRedefine(Cls = y_idx,
                                  OldLabels = 1:nbinsY,
                                  NewLabels = nbinsY:1)
  # print(edges_x)
  # print(table(x_idx))
  
  # Flatten 2D index into 1D: row-major order
  # If y is rows and x is columns (like image/matrix layout)
  flat_idx = (x_idx - 1) * nbinsY + y_idx
  
  # Use tabulate with total number of 2D bins
  counts_flat = tabulate(flat_idx, nbins = nbinsX * nbinsY)
  #print(counts_flat)
  # Reshape to 2D matrix that stores the counted events
  counts_matrix = matrix(counts_flat, nrow = nbinsY, byrow = byrow)
  
  if (isTRUE(names)) {
    colnames(counts_matrix) <- edges_x[1:nbinsX]
    if (isFALSE(redefine)) {
      rownames(counts_matrix) <- edges_y[1:nbinsY]
    } else{
      rownames(counts_matrix) <- rev(edges_y[1:nbinsY])
    }
  }
  
  if (isFALSE(extendOutput)) {
    #default
    return(counts_matrix)
  } else{
    #for tiles
    return(list(
      counts_matrix = counts_matrix,
      x_idx = x_idx,
      y_idx = y_idx
    ))
  }
  
}