#' Local inhomogeneous Spatio-temporal K-functions on a linear network
#'
#' @description
#' The functions \code{localSTLKinhom} and \code{localSTLginhom} implement the
#' inhomogeneous LISTA functions proposed in D'Angelo et al. (2022).
#'
#' @details
#'
#' The homogeneous K-function and pair correlation functions, in
#' D'Angelo et al. (2021), can be obtained easily with \code{localSTLKinhom} and
#' \code{localSTLginhom}, by imputing a lambda vector of constant intensity
#' values, the same for each point.
#'
#' @param x A realisation of a spatio-temporal point processes on a linear network in \code{stlp} format
#' @param lambda values of estimated intensity.
#' @param normalize normalization factor to be considered.
#' @param r values of argument r where K-function will be evaluated. optional.
#' @param t values of argument t where K-function will be evaluated. optional.
#' @param nxy pixel array dimensions. optional.
#' @import spatstat stlnpp
#' @return A list of class \code{lista}.
#' The objects are of class \code{sumstlpp} (Moradi and Mateu, 2020).
#' @export
#'
#' @author Nicoletta D'Angelo
#'
#' @seealso \link{localSTLginhom}, \link[stlnpp]{STLKinhom}, \link[stlnpp]{STLginhom}
#'
#'
#' @examples
#'
#' set.seed(2)
#' df_net <- data.frame(x = runif(25, 0, 0.85), y = runif(25, 0, 0.85), t = runif(25))
#' stlp1 <- stp(df_net, L = chicagonet)
#' lambda <- rep(diff(range(stlp1$df$x)) * diff(range(stlp1$df$y))
#'  * diff(range(stlp1$df$t)) / spatstat.geom::volume(stlp1$L),
#' nrow(stlp1$df))
#'
#' k <- localSTLKinhom(stlp1, lambda = lambda, normalize = TRUE)
#'
#'
#'
#' @references
#' D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' D’Angelo, N., Adelfio, G.,  and Mateu, J. (2022). Local inhomogeneous second-order characteristics for spatio-temporal point processes on linear networks. Stat Papers. https://doi.org/10.1007/s00362-022-01338-4
#'
localSTLKinhom <- function(x, lambda = lambda, normalize = FALSE, r = NULL, t = NULL, nxy = 10)
{
  if (!inherits(x, "stlp"))
    stop("x should be from class stlp")
  
  nx <- nrow(x$df)
  
  if (!is.numeric(lambda)) {
    stop("lambda should be a numeric vector")
  } else {
    if(length(lambda) != nx) {
      stop("lambda should be the same length as the number of points in x")
    }
  } 
  
  if (!is.null(r)){
    if(!is.numeric(r)) {
      stop("r should be a numeric vector")
    } 
  }
  
  if (!is.null(t)){
    if(!is.numeric(t)) {
      stop("t should be a numeric vector")
    } 
  }
  
  if (!is.numeric(nxy)) {
    stop("nxy should be a numeric value")
  } else {
    if(nxy <= 0) {
      stop("nxy should be nxy > 0")
    }
  } 
  
  x$L <- spatstat.linnet::as.linnet(x$L, sparse = FALSE)

  X <- as.stlpp(x)

  Y <- stlnpp::as.lpp.stlpp(X)
  l <- spatstat.geom::domain(Y)
  tleng <- summary(l)$totlength
  n <- spatstat.geom::npoints(X)
  a <- X$time[1]
  b <- X$time[2]
  trange <- b - a
  timev <- X$data$t
  sdist <- spatstat.linnet::pairdist.lpp(Y)
  tdist <- as.matrix(dist(timev))
  toler <- spatstat.linnet::default.linnet.tolerance(l)
  ml <- matrix(1, n, n)
  for (j in 1:n) {
    ml[-j, j] <- spatstat.linnet::countends(l, Y[-j], sdist[-j, j], toler = toler)
  }
  mtplus <- matrix(timev, n, n, byrow = T) + tdist
  mtminus <- matrix(timev, n, n, byrow = T) - tdist
  mtedge <- (mtplus <= b) + (mtminus >= a)
  diag(mtedge) <- 1
  lamden <- outer(lambda, lambda, FUN = "*")
  diag(lamden) <- 1
  edgetl <- mtedge * ml * lamden

  maxs <- 0.7 * max(sdist[!is.infinite(sdist)])
  maxt <- 0.7 * (trange/2)
  if (is.null(r))
    r <- seq((maxs/nxy), maxs, by = (maxs - (maxs/nxy))/(nxy - 1))
  if (is.null(t))
    t <- seq((maxt/nxy), maxt, by = (maxt - (maxt/nxy))/(nxy - 1))

  K_local <- list()
  for(k in 1:n){
    K <- matrix(NA, nrow = nxy, ncol = nxy)
    for (i in 1:length(r)) {
      for (j in 1:length(t)) {
        out <- (sdist <= r[i]) * (tdist <= t[j])
        diag(out) <- 0
        kout <- out /edgetl
        kout <- kout[, k]
        K[i, j] <- sum(kout[!is.na(kout) & !is.infinite(kout)])
      }
    }
    K_local[[k]] <- K
  }

  if (normalize) {
    revrho <- outer(1 / lambda, 1 / lambda, FUN = "*")
    appx <- (n - 1) * (tleng * trange) / (sum(revrho[lower.tri(revrho,diag = FALSE)]) * 2)
    K_local <- lapply(K_local, "*" , appx)

  }
  else {
    K_local <- lapply(K_local, FUN= function(K) (n - 1) * K / (trange * tleng) )
  }
  K_theo <- list()
  for(k in 1:n){
    K_theo[[k]] <- matrix(expand.grid(r, t)[, 1] * expand.grid(r, t)[,2],ncol = nxy)
  }

  Kout <- list()
  for(k in 1:n){
    Kout[[k]] <- list(Kinhom = K_local[[k]], Ktheo = K_theo[[k]], r = r, t = t)
  }

  class(Kout[[k]]) <- c("sumstlpp")
  attr(Kout[[k]], "nxy") <- nxy
  class(Kout) <- "lista"
  return(Kout)
}
