#' Local diagnostics of spatio-temporal point process models
#'
#'
#' @description
#' 
#' This function performs local diagnostics of a model fitted for the
#' first-order intensity of a spatio-temporal point pattern, returning the   points identified as outlying following the diagnostics
#' procedure on individual points of an observed
#' point pattern, as introduced in
#' Adelfio et al. (2020), and applied in D'Angelo et al. (2022) for the linear
#' network case.
#'
#' The points resulting from the local diagnostic procedure provided by this
#' function can be inspected via the \link{plot}, \link{print}, \link{summary},
#' and \link{infl}
#' functions.
#'
#' @details
#' This function performs local diagnostics of a model fitted for the
#' first-order intensity of a spatio-temporal point pattern, by means of  the
#' local spatio-temporal inhomogeneous K-function (Adelfio et al, 2020)
#'  documented by the function
#'  \link[stpp]{KLISTAhat}  of the \code{stpp} package (Gabriel et al, 2013).
#'
#' The function can also perform local diagnostics of a model fitted for the
#' first-order intensity of an  spatio-temporal point pattern on a linear network,
#'  by means of  the
#' local spatio-temporal inhomogeneous K-function on linear networks
#'  (D'Angelo et al, 2021)
#'  documented by the function
#'  \link{localSTLKinhom}.
#'
#' In both cases, it returns the   points identified as outlying following the diagnostics
#' procedure on individual points of an observed
#' point pattern, as introduced in
#' Adelfio et al. (2020), and applied in D'Angelo et al. (2022) for the linear
#' network case.
#' 
#' This function computes discrepancies
#' by means of the \eqn{\chi_i^2} values, obtained following the expression
#' \deqn{
#'   \chi_i^2=\int_L \int_T \Bigg(
#'     \frac{\big(\hat{K}^i_{I}(r,h)- \mathbb{E}[\hat{K}^i(r,h) ]
#'     \big)^2}{\mathbb{E}[\hat{K}^i(r,h) ]}
#'     \Bigg) \text{d}h \text{d}r ,
#'   }
#'     one for each point in the point pattern.
#'     
#' Note that the Euclidean procedure is implemented by the
#' local K-functions of
#'    Adelfio et al. (2020), documented in
#'  \link[stpp]{KLISTAhat} of the \code{stpp} package (Gabriel et al, 2013).
#' The network case uses  the local K-functions on networks (D'Angelo et al., 2021),
#'  documented
#' in  \link{localSTLKinhom}.
#'
#' @param x Either a \code{stp} or a \code{stlp} object
#' @param intensity A vector of intensity values, of the same length as the number
#' of point in \code{x}
#' @param p The percentile to consider as threshold for the outlying points.
#' Default to 0.95.
#'
#' @return A list object of class \code{localdiag}, containing
#' \describe{
#' \item{\code{x}}{The \code{stp} object provided as input}
#' \item{\code{listas}}{The LISTA functions, in a list object}
#' \item{\code{ids}}{The ids of the points identified as outlying}
#' \item{\code{x2}}{A vector with the individual contributions to the Chi-squared statistics,
#'  normalized}
#' \item{\code{p}}{The percentile considered}
#' }
#'
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Giada Adelfio
#'
#' @seealso
#' \link{infl}, \link{plot.localdiag}, \link{print.localdiag},
#' \link{summary.localdiag},
#' \link{globaldiag}
#'
#'
#' @examples
#'
#' set.seed(2)
#' inh <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, 
#'              par = c(.3, 6))
#' 
#' mod1 <- stppm(inh, formula = ~ 1)
#' 
#' resmod1 <- localdiag(inh, mod1$l, p = .9)
#' 
#'
#'
#' @references
#' Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.
#'
#' D’Angelo, N., Adelfio, G.  and Mateu, J. (2022) Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
#' Gabriel, E., Rowlingson, B. S., and Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
#'
localdiag <- function(x, intensity, p = 0.95){
  
  if (!inherits(x, c("stp", "stlp"))) stop("x should be either from class stp or stlp")

  nx <- nrow(x$df)
  
  if (!is.numeric(intensity)) {
    stop("intensity should be a numeric vector")
  } else {
    if(length(intensity) != nx) {
      stop("intensity should be the same length as the number of points in x")
    }
  } 
  
  if (!is.numeric(p)) {
    stop("p should be a numeric value")
  } else {
    if(p > 1 | p < 0) {
      stop("p should be a probability and therefore 0 <= p <= 1")
    }
  } 

  if(any(class(x) == "stp")){
    u <- seq(min(x$df$x, x$df$y), min(x$df$x, x$df$y) / 4, length = 51)
    v <- seq(min(x$df$t), (max(x$df$t) - min(x$df$t)) / 4, length = 51)
    correction_g = "isotropic"
    local_k_1 <- KLISTAhat(xyt = as.stpp(x), dist = u, times = v,
                           correction = correction_g, lambda = intensity)
    theo_local_k <- local_k_1$klistatheo
  } else if(!any(class(x) == "stp")) {
    local_k_1 <- localSTLKinhom(x = x, lambda = intensity)
    theo_local_k <- local_k_1[[1]]$Ktheo
  }


  if(any(class(x) == "stp")) diff_1_local <- array(NA, c(length(u), length(v), nx))

  x2 <- vector(length = nx)

  if(any(class(x) == "stp")){
    for(i in 1:nx){
      diff_1_local[, , i] <- local_k_1$list.KLISTA[, , i] - theo_local_k#
      x2[i] <- sum(diff_1_local[, , i] ^ 2 / (theo_local_k + 0.1))
    }
  } else if(!any(class(x) == "stp")) {
    for(i in 1:nx){
      x2[i] <- sum((local_k_1[[i]]$Kinhom - theo_local_k ^ 2) / theo_local_k)
    }
  }


  normalize <- function(x){(x - min(x)) / ((max(x) - min(x)) + 0.0001) }

  x2 <- normalize(x2)

  id_points <- which(x2 > quantile(x2, p = p))

  listas <- if(any(class(x) == "stp")){
    local_k_1$list.KLISTA
  } else if(!any(class(x) == "stp")) {
    local_k_1
  }

  out <- if(any(class(x) == "stp")){
    list(X = x, listas = listas, ids = id_points, x2 = x2, p = p, 
         dist = local_k_1$dist, times = local_k_1$times)
  } else if(!any(class(x) == "stp")) {
    list(X = x, listas = listas, ids = id_points, x2 = x2, p = p, 
         dist = local_k_1[[1]]$r, times = local_k_1[[1]]$t)
  }
  class(out) <- "localdiag"

  return(out)

}
