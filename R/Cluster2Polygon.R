Cluster2Polygon=function(Data,Cls,ClusterNo,Ratio=1,PlotIt=FALSE){
# polygonzug = Cluster2Polygon(Data,Cls,ClusterNo,Ratio = 1,PlotIt = FALSE)
#
# Constructs a convex or concave geometric boundary for one cluster
# based on the first two columns of Data.
#
# INPUT
# Data        Matrix or data.frame containing the observations.
#             Rows correspond to observations and columns to variables.
#             Data must contain at least two columns. Only the first two columns
#             are used for construction of the geometry and both must be numeric.
#
# Cls
#             Vector of cluster assignments. Its length must be equal to nrow(Data).
#             Missing cluster assignments are ignored.
#
# ClusterNo   Scalar cluster identifier defining which observations are used.
#             Comparison between Cls and ClusterNo is performed after conversion
#             to character.
#
# OPTIONAL
# Ratio       Numeric scalar in [0,1] controlling the geometry.
#             Ratio = 1 returns the convex hull.
#             Ratio < 1 requests a concave hull using sf::st_concave_hull().
#             Smaller values generally allow a more concave boundary.
#             Default: 1.
#
# PlotIt      Logical scalar. If TRUE, the selected observations and the resulting
#             geometry are plotted. Default: FALSE.
#
# OUTPUT
# List V containing exactly four elements:
#        1. Vertices        Numeric matrix with columns:
#                           VertexNo
#                           X
#                           Y
#       2. ClusterNo        Cluster identifier.
#       3. Ratio            Ratio used to construct the geometry.
#       4.                  GeometryType  Geometry type returned for the selected cluster. Possible values
#                           include "POINT", "LINESTRING", and "POLYGON".
#
# Author:  Michael Thrun
  
# DETAILS
# Observations are selected by
#   !is.na(Cls) & as.character(Cls) == as.character(ClusterNo)
# Only the first two columns of Data are used. Rows containing NA, NaN,
# Inf, or -Inf in either of these columns are removed before construction
# of the geometry.
#
# A convex hull is always computed first and is also used to detect
# degenerate clusters:
#   1 hull point  -> POINT
#   2 hull points -> LINESTRING
#   >=3 points    -> POLYGON
# If Ratio = 1, the convex hull is returned.
# If Ratio < 1 and at least three non-collinear hull points are available,
# a concave hull is generated with sf::st_concave_hull(). The package
# 'sf' is therefore required only for non-degenerate concave hulls.
# If Ratio < 1 but the selected cluster is degenerate, no concave hull can
# be constructed. The corresponding POINT or LINESTRING representation is
# returned instead and a warning is issued.
# For polygon geometries generated from the convex hull, the first vertex
# is appended to the end of the path so that the polygon is explicitly
# closed.
  #----------------------------------------------------------------------
  # 1. Input validation
  #----------------------------------------------------------------------
  
  if (is.null(dim(Data)) || ncol(Data) < 2L) {
    stop("Data must contain at least two columns.")
  }
  
  if (length(Cls) != nrow(Data)) {
    stop("length(Cls) must be equal to nrow(Data).")
  }
  
  if (length(ClusterNo) != 1L || is.na(ClusterNo)) {
    stop("ClusterNo must contain exactly one non-missing value.")
  }
  
  if (length(Ratio) != 1L ||
      !is.numeric(Ratio) ||
      !is.finite(Ratio) ||
      Ratio < 0 ||
      Ratio > 1) {
    stop("Ratio must be exactly one finite numeric value in [0, 1].")
  }
  #----------------------------------------------------------------------
  # 2. Select observations belonging to ClusterNo
  #----------------------------------------------------------------------
  
  keep = !is.na(Cls) & as.character(Cls) == as.character(ClusterNo)
  
  XY = Data[keep, 1:2, drop = FALSE]
  
  if (!is.numeric(XY[, 1]) || !is.numeric(XY[, 2])) {
    stop("The first two columns of Data must be numeric.")
  }
  
  #----------------------------------------------------------------------
  # 3. Remove non-finite and duplicate coordinate pairs
  #----------------------------------------------------------------------
  
  # chull() akzeptiert keine NA-, NaN- oder Inf-Werte.
  XY <- XY[
    is.finite(XY[,1]) & is.finite(XY[,2]),
    ,
    drop = FALSE
  ]
  #unique kann nun in R die punkte unqique machen!
  XY=unique(XY)
  
  rownames(XY) = NULL
  
  if (nrow(XY) == 0L) {
    stop(
      "Cluster ",
      as.character(ClusterNo),
      " contains no finite observations in the first two columns of Data."
    )
  }
  #------------------------------------------------------------
  # 4. Compute convex reference hull
  #
  # Sie dient zugleich dazu, degenerierte Cluster zu erkennen:
  #   1 Hüllpunkt  -> POINT
  #   2 Hüllpunkte -> LINESTRING / kollinear
  #   >= 3         -> POLYGON
  #------------------------------------------------------------∂
  
  if (nrow(XY) <= 2L) {
    
    convex_index = seq_len(
      nrow(XY)
    )
    
  } else {
    
    convex_index = grDevices::chull(
      XY[, 1L],
      XY[, 2L]
    )
  }
  #----------------------------------------------------------------------
  # 5. Convex hull or degenerate geometry
  #----------------------------------------------------------------------
  is_degenerate = length(convex_index) < 3L
  if(Ratio==1 ||is_degenerate){ #konvex 
    geometry_type = if (length(convex_index) == 1L) {
      "POINT"
    } else if (length(convex_index) == 2L) {
      "LINESTRING"
    } else {
      "POLYGON"
    }
    # Only actual polygons are explicitly closed.
    if (geometry_type == "POLYGON") {
      
      path_index = c(
        convex_index,
        convex_index[1L]
      )
      
    } else {
      
      path_index = convex_index
    }
    
    path = XY[
      path_index,
      1:2,
      drop = FALSE
    ]
    
    
    if (Ratio < 1 && is_degenerate) {
      
      warning(
        "Cluster ",
        as.character(ClusterNo),
        " has fewer than three non-collinear hull points. ",
        "A concave polygon hull is not defined; ",
        geometry_type,
        " is returned instead."
      )
    }
    
    
    #----------------------------------------------------------------------
    # 6. Concave hull
    #----------------------------------------------------------------------
    
  } else {
    
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop(
        "For Ratio < 1 the package 'sf' is required."
      )
    }
    
    
    # Construct exactly one MULTIPOINT geometry.
    cluster_geometry = sf::st_sfc(
      sf::st_multipoint(
        as.matrix(XY)
      )
    )
    
    
    hull = sf::st_concave_hull(
      cluster_geometry,
      ratio = Ratio,
      allow_holes = FALSE
    )
    
    
    coordinates =
      sf::st_coordinates(hull)
    
    
    path = cbind(
      X = coordinates[, "X"],
      Y = coordinates[, "Y"]
    )
    
    
    geometry_type = as.character(
      sf::st_geometry_type(hull)[1L]
    )
  }
  
  
  #----------------------------------------------------------------------
  # 7. Construct vertex matrix
  #----------------------------------------------------------------------
  
  vertex_matrix = cbind(
    
    VertexNo = seq_len(
      nrow(path)
    ),
    
    X = as.numeric(
      path[, 1L]
    ),
    
    Y = as.numeric(
      path[, 2L]
    )
  )
  
  rownames(vertex_matrix) = NULL
  
  
  #----------------------------------------------------------------------
  # 8. Construct requested four-element output list
  #----------------------------------------------------------------------
  
  polygonzug = list(
    
    Vertices = vertex_matrix,
    
    ClusterNo = ClusterNo,
    
    Ratio = Ratio,
    
    GeometryType = geometry_type
  )
  
  
  #----------------------------------------------------------------------
  # 9. Optional plot
  #----------------------------------------------------------------------
  
  if (isTRUE(PlotIt)) {
    
    axis_names = colnames(Data)
    
    if (is.null(axis_names)) {
      axis_names = c(
        "X",
        "Y"
      )
    }
    
    
    graphics::plot(
      XY[, 1L],
      XY[, 2L],
      xlab = axis_names[1L],
      ylab = axis_names[2L],
      pch = 16,
      main = paste0(
        "Cluster ",
        as.character(ClusterNo),
        ", Ratio = ",
        Ratio
      )
    )
    
    
    if (nrow(vertex_matrix) >= 2L) {
      
      graphics::lines(
        vertex_matrix[, "X"],
        vertex_matrix[, "Y"],
        col = "red",
        lwd = 2
      )
      
    } else {
      
      graphics::points(
        vertex_matrix[, "X"],
        vertex_matrix[, "Y"],
        col = "red",
        pch = 16,
        cex = 1.5
      )
    }
  }
  
  
  return(polygonzug)
}