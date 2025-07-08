#' Plot of the global diagnostics of a spatio-temporal point process first-order intensity
#'
#' @description
#' This function performs global diagnostics of a model fitted for the
#' first-order intensity of a  spatio-temporal point pattern, by returning
#' the plots of the inhomogeneous K-function weighted by the
#' provided intensity to diagnose, its theoretical value,
#' and their difference.
#'
#'
#' @param x A \code{globaldiag} object
#' @param samescale Logical value. It indicates whether to plot the observed
#' and the theoretical K-function in the same or
#' different scale. Default to \code{TRUE}.
#' @param ... additional unused argument
#'
#'
#' @return It plots three panels: the observed K-function, as returned by \link[stlnpp]{STLKinhom};
#' the theoretical one; their difference. The function also prints the sum of
#'  squared differences between the observed and theoretical
#' K-function on the console.
#'
#' @export
#'
#' @author Nicoletta D'Angelo 
#'
#' @seealso
#' \link{globaldiag}, \link{print.globaldiag},
#' \link{summary.globaldiag}
#'
#'
#'
#'
#' @examples
#'
#' set.seed(2)
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'                par = c(.3, 6))
#' 
#' mod1 <- stppm(inh, formula = ~ 1)
#' mod2 <- stppm(inh, formula = ~ x)
#' 
#' g1 <- globaldiag(inh, mod1$l)
#' g2 <- globaldiag(inh, mod2$l)
#' 
#' plot(g1)
#' plot(g2)
#'
#' @references
#' Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.
#'
#' D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
#' Gabriel, E., and Diggle, P. J. (2009). Second‐order analysis of inhomogeneous spatio‐temporal point process data. Statistica Neerlandica, 63(1), 43-51.
#'
#' Gabriel, E., Rowlingson, B. S., & Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
#'
#' Moradi M, Cronie O, and Mateu J (2020). stlnpp: Spatio-temporal analysis of point patterns on linear networks.
#'
#' Moradi, M. M., and Mateu, J. (2020). First-and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics, 29(3), 432-443.
#'
#'
plot.globaldiag <- function(x, samescale = TRUE, ...){
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  lims <- if(samescale){
    range(c(as.numeric(x$est), as.numeric(x$theo), as.numeric(x$diffK)))
  } else {
    NULL
  }
  
  par(mfrow = c(1, 3))
  
  fields::image.plot(x$dist, x$times, x$est, main = "Estimated", xlab = "r",
                     ylab = "h",
                     col = grDevices::hcl.colors(12, "YlOrRd", rev = TRUE),
                     zlim = lims,
                     legend.mar =  12, axes = FALSE)
  axis(1, at = seq(0, 1, l = length(x$dist)), labels = round(x$dist, 3))
  axis(2, at = seq(0, 1, l = length(x$times)), labels = round(x$times, 3))
  box()
  fields::image.plot(x$dist, x$times, x$theo, main = "Theoretical", xlab = "r",
                     ylab = "h",
                     col = grDevices::hcl.colors(12, "YlOrRd", rev = TRUE),
                     zlim = lims,
                     legend.mar =  12, axes = FALSE)
  axis(1, at = seq(0, 1, l = length(x$dist)), labels = round(x$dist, 3))
  axis(2, at = seq(0, 1, l = length(x$times)), labels = round(x$times, 3))
  box()
  
  fields::image.plot(x$dist, x$times, x$diffK, main = "Difference", xlab = "r",
                     ylab = "h",
                     col = grDevices::hcl.colors(12, "YlOrRd", rev = TRUE),
                     zlim = lims,
                     legend.mar =  15, axes = FALSE)
  axis(1, at = seq(0, 1, l = length(x$dist)), labels = round(x$dist, 3))
  axis(2, at = seq(0, 1, l = length(x$times)), labels = round(x$times, 3))
  box()
  
}
