#' isomap function
#'
#' @param sf TODO
#' @param values TODO
#' @param breaks TODO
#' @param x.grid TODO
#' @param y.grid TODO
#' @param pr.grid TODO
#' @param cells TODO
#' @param n TODO
#' @param name TODO
#' @param divcol TODO
#' @param placelegend TODO
#' @param space TODO
#' @param showlegend TODO
#' @param ... TODO
#'
#' @return TODO
#' @export

isomap <- function(sf, values, breaks, x.grid, y.grid, pr.grid, cells, n = 9, name = "PuBu", divcol = FALSE, placelegend = "right", space = 8, showlegend = TRUE, ...){

  class <- classInt::classIntervals(values,
                                    style = "fixed",
                                    fixedBreaks = breaks)
  pal <- RColorBrewer::brewer.pal(n, name)
  if (divcol) {pal <- rev(pal)}
  col.pal <- colorRampPalette(pal)((length(breaks)-1))
  which.col <- classInt::findColours(class, col.pal)

  Show <- rep(NA, length(x.grid))
  Show[pr.grid > 0] <- values
  Show <- matrix(Show, cells[1], cells[2])

  op <- par(
    oma = c(0, 0, 0, space), # Room for title and legend
    mfrow = c(1, 1) # Add more plots if you want
  )

  #   image(unique(x.grid), unique(y.grid), Show,
  #         col = col.pal,
  #         breaks = breaks,
  #         xlab = "", ylab = "", axes = FALSE)
  #   plot(sf, add = TRUE, ...)

  plot(sf)
  image(unique(x.grid), unique(y.grid), Show,
        col = col.pal,
        breaks = breaks,
        xlab = "", ylab = "", axes = FALSE, add = TRUE)
  plot(sf, col = NA, add = TRUE, ...)

  par(op) # Leave the last plot

  op <- par(usr = c(0, 1, 0, 1), # Reset the coordinates
            xpd = NA) # Allow plotting outside the plot region

  if (showlegend) {
    legend(placelegend, bty = "n",
           fill = rev(attr(which.col, "palette")),
           legend = rev(names(attr(which.col, "table"))))
  }

}