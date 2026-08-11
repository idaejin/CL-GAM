#' atp_setup2
#'
#' @param sf TODO
#' @param coords TODO
#' @param y TODO
#' @param e TODO
#' @param Cweights TODO
#' @param show TODO
#'
#' @return TODO
#' @export

atp_setup2 <- function(sf, coords, y, e, Cweights = FALSE, show = FALSE){
  # Obtain shapefile's polygons
  sf.poly <- sp::SpatialPolygons(sf@polygons)

  # Total number of polygons
  n.poly <- length(sf.poly)

  # Set grid (centroids of each cell)
  x.grid <- coords[,1]
  y.grid <- coords[,2]

  # Determine the grid points for estimation that fall inside map
  pr.grid <- 0 * x.grid

  for (k in 1:n.poly) {
    # Points that fall in polygon k
    pp <- over(sp::SpatialPoints(coords), sf.poly[k]) == 1
    pr.grid[pp] <- k
  }

  selp <- pr.grid > 0
  x.pred <- x.grid[selp]
  y.pred <- y.grid[selp]
  pr.pred <- pr.grid[selp]

  ### Obtain spatial composition matrix
  C <- spam::spam(list(i = pr.pred, j = 1:(length(pr.pred)), x = rep(1, length(pr.pred))))
  C <- as.matrix(C)

  ### Obtain naive estimates for the grid points
  if (any(rowSums(C) == 0)) stop("Increase the number of cells")
  C.star <- (1/rowSums(C)) * C
  y.naive <- c(crossprod(C.star, y))
  e.naive <- c(crossprod(C.star, e))

  ### Output (first option)
  list1 <- list(y.naive = y.naive, e.naive = e.naive,
                C = C, x.pred = x.pred, y.pred = y.pred,
                pr.pred = pr.pred, x.grid = x.grid,
                y.grid = y.grid, pr.grid = pr.grid)

  if (Cweights) {
    # Convert to 'gpc.poly' format
    map.gpc <- as(sf, 'gpc.poly')
    # Obtain the neighbours of each polygon (in matrix format)
    adj <- rgeos::gTouches(sf, byid = TRUE)
    # Setup to construct the cells of each prediction point
    dx <- 0.5 * diff(xdiv)[1]
    dy <- 0.5 * diff(ydiv)[1]
    # *Note: There exist a very few discrepancy between the
    #        elements of diff(xdiv) (and diff(y div) also)!
    #        So, I pick up the first element of each vector
    ddx <- c(-1, 1, 1, -1) * dx
    ddy <- c(-1, -1, 1, 1) * dy
    # *Note: ddx and ddy will help me to construct the vertices
    #        of each cell
    # Set new composition matrix
    C2 <- C
    # *Note: I only modify some elements of C, using the
    #        following loop ...

    for (k in 1:n.poly) {
      # For each polygon, select the points that fall inside it
      xx <- x.pred[pr.pred == k]
      yy <- y.pred[pr.pred == k]
      # How many points?
      npoints <- length(xx)
      # Pick up the selected polygon
      P <- map.gpc[[k]]
      # Pick up the correponding neighbours
      neigh <- which(adj[k,] == TRUE)
      nneigh <- length(neigh)
      list.neigh <- map.gpc[neigh]

      for (j in 1:npoints) {
        # For each point, construct its cell and calculate
        # the corresponding area
        xr <- xx[j] + ddx
        yr <- yy[j] + ddy
        R <- as(cbind(xr,  yr), 'gpc.poly')
        ra <- rgeos::area.poly(R)
        # Calculate intersection between P and R
        S <- rgeos::intersect(P, R)
        sa <- rgeos::area.poly(S)
        # Find the column position of this point in
        # the composition matrix C
        colpos <- which(x.pred == xx[j] & y.pred == yy[j])
        accum <- 0
        for (i in 1:nneigh){
          # For each neighbour, obtain the area of its
          # intersection with R
          aux <- rgeos::intersect(list.neigh[[i]], R)
          #plot(aux, poly.args = list(col ='green'), add = TRUE)
          int.neigh <- rgeos::area.poly(aux)
          int.area <- int.neigh/ra
          # Put this area in C2 (if appropiate)
          if (sa < 0.999*ra && int.neigh != 0){
            C2[neigh[i], colpos] <- int.area
          }
          # Accumulate these areas
          accum <- accum + int.area
        }
        # Put the remaining area in the position of the cell
        C2[k, colpos] <- 1 - accum
      }
    }

    # Normalize the columns that do not sum up 1, as follows:
    which.norm <- which(colSums(C2) != 1)
    for (j in which.norm) {
      # For the columns given in which.norm, normalize them as:
      t <- sum(C2[, j])
      max.pos <- which(C2[, j] == max(C2[, j]))
      C2[max.pos, j] <- C2[max.pos, j] + (1-t)
    }
    # Obtain new naive estimates at the prediction grid points
    # inside the map
    C.star.w <- (1/rowSums(C2)) * C2
    y.naive.w <- c(crossprod(C.star.w, y))
    e.naive.w <- c(crossprod(C.star.w, e))

    ### Output (second option)
    list2 <- list(y.naive.w = y.naive.w, e.naive.w = e.naive.w, C.w = C2)
    list3 <- c(list1, list2)

  }

  ### Show grid points for estimation
  if (show == TRUE) {
    x11(title = "Fine grid for estimation (blue points)")
    plot(sf)
    points(cbind(x.pred, y.pred), col = "blue", cex = 0.1,
           pch = 19)
    points(coordinates(sf), col = "red", cex = 0.5,
           pch = 19)
  }

  ifelse(Cweights == FALSE, return(list1), return(list3))
}