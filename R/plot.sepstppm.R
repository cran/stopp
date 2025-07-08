#' Plot of the fitted intensity of a separable spatio-temporal Poisson model
#'
#' The function plots the fitted intensity, displayed both in space and in space and time.
#'
#' @param x An object of class \code{sepstppm}
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
#' @param par Default to \code{TRUE}. If \code{par=FALSE},
#' the user is asked for input, before a new figure is drawn.
#' @param sig Smoothing bandwidth for spatial representation
#' @param ... additional unused argument
#'
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#'
#' @examples
#'
#' crimesub <- stpm(valenciacrimes$df[1:100, ],
#'            names = colnames(valenciacrimes$df)[-c(1:3)])
#' 
#' mod1 <- sepstppm(crimesub, spaceformula = ~x ,
#'                   timeformula = ~ day)
#'                   
#' plot(mod1)
#'
plot.sepstppm <- function(x,
                          scaler = c("silverman", "IQR", "sd", "var"),
                          do.points = TRUE,
                          print.bw = FALSE,
                          zap = 0.00001,
                          par = TRUE,
                          sig = NULL,
                          ...){
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  mark_int <- x$l
  
  ppx_int <- suppressWarnings(spatstat.geom::ppp(x$x$df$x, x$x$df$y, marks = mark_int,
                                                 spatstat.geom::owin(range(x$x$df$x),
                                                                     range(x$x$df$y))))
  
  sig <- sparr::OS(unmark(ppx_int), scaler = scaler)
  
  if(par == T){
    par(mfrow = c(1, 2))
    par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, 1))
    plot(spatstat.explore::Smooth(ppx_int, sigma = sig), zap = zap,
         col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                             range = range(mark_int)),
                    "stuff")$outputs,
         main = c("Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
    
    par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))
    
    plot3D::scatter3D(x$x$df$x, x$x$df$y, x$x$df$t,
                      theta = - 45, phi = 20,
                      col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                          range = range(mark_int)),
                                 "stuff")$outputs,
                      ticktype = "detailed", pch = 20,
                      colvar = mark_int,
                      xlab="x",ylab="y",zlab="t",
                      main = c("Intensity in space-time \n Pointwise computation"))
    par(mar = c(5, 4, 4, 2) + 0.1)
  } else {
    par(mfrow = c(1, 1))
    par(ask = TRUE)
    plot(spatstat.explore::Smooth(ppx_int, sigma = sig), zap = zap,
         col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                             range = range(mark_int)),
                    "stuff")$outputs,
         main = c("Intensity in space \n Density Kernel Smoothing"))
    if(do.points == T){plot(spatstat.geom::unmark(ppx_int), add = T)}
    
    par(ask = TRUE)
    par(mar = c(5, 4, 4, 2) + 0.1 - c(4, 1 , 1, -0.5))
    plot3D::scatter3D(x$x$df$x, x$x$df$y, x$x$df$t,
                      theta = - 45, phi = 20,
                      col = attr(spatstat.geom::colourmap(grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                                                          range = range(mark_int)),
                                 "stuff")$outputs,
                      ticktype = "detailed", pch = 20,
                      colvar = mark_int,
                      xlab="x",ylab="y",zlab="t",
                      main = c("Intensity in space-time \n Pointwise computation"))
    par(mar = c(5, 4, 4, 2) + 0.1)
    par(ask = FALSE)
  }
  if(print.bw == T){print(sig)}
}

