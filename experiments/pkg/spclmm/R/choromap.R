#' choromap function
#'
#' @param sf SpatialPolygons*, sf object, or sfc geometry
#' @param values numeric values aligned with polygons
#' @param breaks break points for classIntervals
#' @param n palette size
#' @param name RColorBrewer palette name
#' @param divcol reverse palette if TRUE
#' @param placelegend legend position
#' @param title plot title
#' @param space oma right margin for legend
#' @param showlegend draw legend
#' @param ... passed to plot / plot.sfc
#'
#' @return invisible NULL
#' @export
choromap <- function(sf, values, breaks, n = 9,
                     name = "PuBu", divcol = FALSE,
                     placelegend = "right",
                     title = "", space = 8,
                     showlegend = TRUE, ...){

  class <- classInt::classIntervals(values,
                                    style = "fixed",
                                    fixedBreaks = breaks)
  pal <- RColorBrewer::brewer.pal(n, name)
  if (divcol) {pal <- rev(pal)}
  col.pal <- colorRampPalette(pal)((length(breaks)-1))
  which.col <- classInt::findColours(class, col.pal)

  # Prefer sf geometry (sp::plot can fail under modern sf/sp stacks)
  geom <- NULL
  if (inherits(sf, "sf")) {
    geom <- sf::st_geometry(sf)
  } else if (inherits(sf, "sfc")) {
    geom <- sf
  } else if (inherits(sf, "Spatial") && requireNamespace("sf", quietly = TRUE)) {
    geom <- sf::st_geometry(sf::st_as_sf(sf))
  }

  op <- par(
    oma = c(0, 0, 0, space),
    mfrow = c(1, 1)
  )
  if (!is.null(geom)) {
    plot(geom, col = which.col, ...)
  } else {
    # legacy fallback
    sp::plot(sf, col = which.col, ...)
  }
  title(main = paste(title))

  par(op)

  op <- par(usr = c(0, 1, 0, 1), xpd = NA)

  if (showlegend) {
    legend(placelegend, bty = "n",
           fill = rev(attr(which.col, "palette")),
           legend = rev(names(attr(which.col, "table"))))
  }
  invisible(NULL)
}
