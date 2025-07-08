#' Plot a \code{stlp} object
#'
#' This function plots the point pattern on a linear network
#'  stored in the \code{stlp} object given in input,
#' in a three panel plot representing the \link[plot3D]{plot3D} of the coordinates, and the
#' marginal spatial and temporal coordinates.
#'
#' @param x An object of class \code{stp}
#' @param tcum If \code{TRUE} (default option),
#' the temporal point pattern is displayed cumulatively.
#' A barplot is automatically plotted if there are repeated
#' counts (typically with discrete times).
#' @param marg Default to \code{TRUE}. If \code{FALSE}, only the
#' spatio-temporal point pattern is plotted.
#' @param col The color of the points. Default to \code{"black"}
#' @param cols The color of the linear network. Default to \code{"grey"}
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso
#' \link{stp}, \link{summary.stlp}, \link{print.stlp}
#'
#' @examples
#'
#'
#' set.seed(2)
#' df_net <- data.frame(cbind(runif(100, 0, 0.85), runif(100, 0, 0.85), runif(100)))
#' 
#' stlp1 <- stp(df_net, L = chicagonet)
#' 
#' plot(stlp1)
#'
plot.stlp <- function(x, tcum = TRUE, marg = TRUE, col = 1, cols = "grey", ...)
{
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if(marg == TRUE) par(mfrow = c(1, 3)) else par(mfrow = c(1, 1))

  plot3D::scatter3D(x$df$x, x$df$y, x$df$t,
                    theta = - 45, phi = 20,  pch = 20,
                    ticktype = "detailed",
                    col = col,
                    xlab = "x", ylab = "y", zlab = "t",
                    main = c("Spatio-temporal point pattern"))

  if(marg == TRUE) {
    plot(x$L, col = cols, main = "Spatial point pattern")
    title(xlab = "x", ylab = "y")
    axis(1)
    axis(2)
    box()
    points(x$df$x, x$df$y, pch = 20, col = col)
    if(any(table(x$df$t) > 1)){
      barplot(table(x$df$t), main = "Temporal point pattern", xlab = "t")
    } else {
      if(tcum == TRUE){
        plot(sort(x$df$t), sort(cumsum(x$df$t)), main = "Temporal point pattern",
             xlab = "t", ylab = "", yaxt = "n", type = "l", pch = 20, col = col)
      }else{
        plot(x$df$t, rep(1, length(x$df$t)), main = "Temporal point pattern",
             xlab = "t", ylab = "", yaxt = "n", pch = 20, col = col)
      }

    }
  }
}
