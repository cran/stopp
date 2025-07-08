#' Plot of the fitted intensity of a LGCP model
#'
#' The function plots the fitted intensity, displayed both in space and in space and time.
#'  In the case of local covariance parameters, the function returns the
#'  mean of the random intensity, displayed both in space and in space and time.
#'
#' @param x An object of class \code{stlgcppm}
#' @param scaler Optional. Controls the value for a scalar representation of the
#'  spatial scale of the data.
#'  Either a character string, \code{"silverman"} (default), \code{"IQR"},
#'   \code{"sd"}, or \code{"var"};
#'  or positive numeric value(s). See \link[sparr]{OS}.
#' @param do.points Add points to plot
#' @param print.bw It prints the estimated oversmoothing (\link[sparr]{OS}) bandwidth selector
#' @param zap Noise threshold factor (default to 0.00001). A numerical value greater than or equal to 1.
#'  If the range of pixel values is less than \code{zap * .Machine$double.eps},
#'   the image will be treated as constant. This avoids displaying images which
#'   should be constant but contain small numerical errors.
#' @param par Default to \code{TRUE}.
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{stlgcppm}, \link{print.stlgcppm}, \link{summary.stlgcppm},
#' \link{localsummary}, \link{localplot}
#'
#'
#'
#' @examples
#' \donttest{
#' catsub <- stp(greececatalog$df[1:200, ])
#' 
#' lgcp_loc <- stlgcppm(catsub, formula = ~ x, first = "local")
#'
#' plot(lgcp_loc)
#' 
#'}
#'
#'
#'
#' @references
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
#' Davies, T.M. and Hazelton, M.L. (2010), Adaptive kernel estimation of spatial relative risk, Statistics in Medicine, 29(23) 2423-2437.
#'
#' Siino, M., Adelfio, G., and Mateu, J. (2018). Joint second-order parameter estimation for spatio-temporal log-Gaussian Cox processes. Stochastic environmental research and risk assessment, 32(12), 3525-3539.
#'
#' Terrell, G.R. (1990). The maximal smoothing principle in density estimation, Journal of the American Statistical Association, 85, 470-477.
#'
#'
#'
plot.stlgcppm <- function(x,
                          scaler = c("silverman", "IQR", "sd", "var"),
                          do.points = TRUE,
                          print.bw = FALSE,
                          zap = 0.00001,
                          par = TRUE,
                          ...){
  
  if(inherits(x$CovCoefs, "numeric") & inherits(x$IntCoefs, "numeric") & length(x$IntCoefs) == 1){
    stop("Constant intensity, no plot to show")
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  mark_int <- x$l
  mark_all <- x$mu

  ppx_int <- suppressWarnings(spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = mark_int,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y))))
  ppx_all <- suppressWarnings(spatstat.geom::ppp(x$X$df$x, x$X$df$y, marks = mark_all,
                                window = spatstat.geom::owin(range(x$X$df$x), range(x$X$df$y))))
  sig <- sparr::OS(unmark(ppx_int), scaler = scaler)

  if(inherits(x$CovCoefs, "numeric")){
    if(par == T){
      par(mfrow = c(1, 2))
      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, 1))
      plot(spatstat.explore::Smooth(ppx_int, sigma = sig), zap = zap,
           col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                               range = range(mark_int)),
                      "stuff")$outputs,
           main = c("First-order Intensity in space \n Density Kernel Smoothing"))
      if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}

      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))

      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(mark_int)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 20,
                        colvar = mark_int,
                        xlab="x",ylab="y",zlab="t",
                        main = c("First-order Intensity in space-time \n Pointwise computation"))
      par(mar = c(5, 4, 4, 2) + 0.1)
    } else {
      par(mfrow = c(1, 1))
      par(ask = TRUE)
      plot(spatstat.explore::Smooth(ppx_int, sigma = sig), zap = zap,
           col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                               range = range(mark_int)),
                      "stuff")$outputs,
           main = c("First-order Intensity in space \n Density Kernel Smoothing"))
      if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
      par(ask = TRUE)
      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(mark_int)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 20,
                        colvar = mark_int,
                        xlab="x",ylab="y",zlab="t",
                        main = c("First-order Intensity in space-time \n Pointwise computation"))
      par(mar = c(5, 4, 4, 2) + 0.1)
      par(ask = FALSE)
    }

  } else {

    if(inherits(x$IntCoefs, "numeric") & length(x$IntCoefs) == 1){
      if(par == T){
        par(mfrow = c(1, 2))
      } else {
        par(mfrow = c(1, 1))
        par(ask = TRUE)
      }

      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, 1))
      plot(spatstat.explore::Smooth(ppx_all, sigma = sig), zap = zap,
           col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                               range = range(mark_all)),
                      "stuff")$outputs,
           main = c("Mean Function of Lambda in space \n Density Kernel Smoothing"))
      if(do.points == T){plot(spatstat.geom::unmark(ppx_all), add = T)}

      if(par != T){

        par(ask = TRUE)
      }

      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(mark_all)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 20,
                        colvar = mark_all,
                        xlab="x",ylab="y",zlab="t",
                        main = c("Mean Function of Lambda in space-time \n Pointwise computation"))
      par(mar = c(5, 4, 4, 2) + 0.1)
      if(par != T){
        par(ask = FALSE)
      }

    } else {

      if(par == T){
        par(mfrow = c(2, 2))
      } else {
        par(mfrow = c(1, 1))
        par(ask = TRUE)
      }

      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, 1))
      plot(spatstat.explore::Smooth(ppx_int, sigma = sig), zap = zap,
           col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                               range = range(mark_int)),
                      "stuff")$outputs,
           main = c("First-order Intensity in space \n Density Kernel Smoothing"))
      if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}

      if(par != T){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(mark_int)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 20,
                        colvar = mark_int,
                        xlab="x",ylab="y",zlab="t",
                        main = c("First-order Intensity in space-time \n Pointwise computation"))

      if(par != T){
        par(ask = TRUE)
      }

      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, 1))
      plot(spatstat.explore::Smooth(ppx_all, sigma = sig), zap = zap,
           col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                               range = range(mark_all)),
                      "stuff")$outputs,
           main = c("Mean Function of Lambda in space \n Density Kernel Smoothing"))
      if(do.points == T){plot(spatstat.geom::unmark(ppx_all), add = T)}

      if(par != T){
        par(ask = TRUE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))
      plot3D::scatter3D(x$X$df$x, x$X$df$y, x$X$df$t,
                        theta = - 45, phi = 20,
                        col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                            range = range(mark_all)),
                                   "stuff")$outputs,
                        ticktype = "detailed", pch = 20,
                        colvar = mark_all,
                        xlab="x",ylab="y",zlab="t",
                        main = c("Mean Function of Lambda in space-time \n Pointwise computation"))
      if(par != T){

        par(ask = FALSE)
      }
      par(mar = c(5, 4, 4, 2) + 0.1)

    }}

  if(print.bw == T){print(sig)}

}
