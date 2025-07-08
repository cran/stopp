#' Test of local structure for spatio-temporal point processes
#'
#' @description
#' This function performs the permutation test of the local structure for
#' spatio-temporal point pattern data, proposed in Siino et al. (2018), as well as for
#' spatio-temporal point pattern data occurring on the same linear network,
#'  following
#' D'Angelo et al. (2021).
#'
#' @param X Background spatio-temporal point pattern. Usually, the most clustered
#' between \code{X} and \code{Z}. Must be either a \code{stp} or \code{stlp} object.
#' @param Z Other spatio-temporal point pattern. Must also be of the same class as \code{X}.
#' @param method Character string indicating which version of LISTA function to use:
#' either
#'  \code{"K"} or \code{"g"}.
#' If \code{"K"}, the local spatio-temporal K-function is used to run the test.
#' If \code{"g"}, the local spatio-temporal pair correlation function is used.
#' @param k Number of permutations
#' @param alpha Significance level
#' @param verbose If TRUE (default) the progress of the test is printed
#'
#'
#' @details
#' The test detects local differences between  \eqn{\textbf{x}} and \eqn{\textbf{z}}
#' occurring on the same space-time region.
#'
#'  The test ends providing a vector \eqn{p} of  \eqn{p}- values, one for each point
#'  in \eqn{\textbf{x}}.
#'
#'  If the test is  performed for spatio-temporal point patterns as in
#'    Siino et al. (2018), that is, on an object of class \code{stp}, the LISTA
#'    functions \eqn{\hat{L}^{(i)}} employed are the local functions of
#'    Adelfio et al. (2020), documented in
#'  \link[stpp]{KLISTAhat} and \link[stpp]{LISTAhat} of the \code{stpp} package (Gabriel et al, 2013).
#'
#'  If the function is applied to a \code{stlp} object, that is, on two spatio-temporal
#'  point patterns observed on the same linear network \code{L}, the LISTA function
#' \eqn{\hat{L}^{(i)}} used are the ones proposed in D'Angelo et al. (2021), documented
#' in  \link{localSTLKinhom} and \link{localSTLginhom}.
#'
#' Details on the performance of the test are found in Siino et al. (2018) and
#' D'Angelo et al. (2021), for Euclidean and network spaces, respectively.
#'
#'
#' @return A list of class \code{localtest}, containing
#' \describe{
#' \item{\code{p}}{A vector of p-values, one for each of the points in \code{X}}
#' \item{\code{X}}{The background spatio-temporal point pattern given in input}
#' \item{\code{Z}}{The alternative spatio-temporal point pattern given in input}
#' \item{\code{alpha}}{The threshold given in input}
#' \item{\code{Xsig}}{A \code{stp} object storing the resulting significant points}
#' \item{\code{Xnosig}}{A \code{stp} object storing the resulting non-significant points}
#' \item{\code{id}}{The ids of the resulting significant points}
#' }
#'
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Marianna Siino
#'
#' @seealso
#' \link{print.localtest}, \link{summary.localtest}, \link{plot.localtest}
#'
#'
#' @examples
#'
#'
#'
#'\donttest{
#'
#' set.seed(2)
#' X <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)},
#'             par = c(.005, 5))
#' Z <- rstpp(lambda = 30)
#' 
#' test <- localtest(X, Z, method = "K", k = 3)
#'
#'}
#'
#'
#' @references
#' Adelfio, G., Siino, M., Mateu, J., and Rodríguez-Cortés, F. J. (2020). Some properties of local weighted second-order statistics for spatio-temporal point processes. Stochastic Environmental Research and Risk Assessment, 34(1), 149-168.
#'
#' D’Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.
#'
#' Gabriel, E., Rowlingson, B. S., and Diggle, P. J. (2013). stpp: An R Package for Plotting, Simulating and Analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1–29. https://doi.org/10.18637/jss.v053.i02
#'
#' Siino, M., Rodríguez‐Cortés, F. J., Mateu, J. ,and Adelfio, G. (2018). Testing for local structure in spatiotemporal point pattern data. Environmetrics, 29(5-6), e2463.
#'
localtest <- function(X, Z, method = c("K", "g"), k, alpha = 0.05, verbose = TRUE){
  
  if (!inherits(X, c("stp", "stlp"))) stop("X should be either from class stp or stlp")
  
  if (!is.numeric(k)) {
    stop("k should be a numeric value")
  } else {
    if(k <= 2) {
      stop("k should be k >= 3")
    }
  } 

  if (!is.numeric(alpha)) {
    stop("alpha should be a numeric value")
  } else {
    if(alpha > 1 | alpha < 0) {
      stop("alpha should be a probability and therefore 0 <= alpha <= 1")
    }
  } 
  
  method <- match.arg(method)
  
  nX <- nrow(X$df)
  nZ <- nrow(Z$df)
  perm <- c(rep(1, nZ), rep(2, nX - 1))
  p <- double(nX)

  if ((inherits(X,"stlp")|inherits(Z,"stlp"))) L0 <- X$L

  res_test <- double(nX)
  res_test_d <- matrix(NA, k, nX)

  if ((inherits(X,"stp")|inherits(Z,"stp"))){

    X_lista <- if (method == "K"){
      stpp::KLISTAhat(as.stpp(X))$list.KLISTA
    } else if (method == "g"){
      stpp::LISTAhat(as.stpp(X))$list.LISTA
    } else {stop(" 'method' argument must be either \"g\" or \"K\ ")}

  } else if ((inherits(X,"stlp")|inherits(Z,"stlp"))){
    X_lista <- if (method == "K"){
      STLKinhom_i(as.stlpp(X), lambda = rep(nX / (spatstat.geom::volume(L0) * (range(X$df$t)[2] - range(X$df$t)[1])), nX))
    } else if (method == "g"){
      STLginhom_i(as.stlpp(X), lambda = rep(nX / (spatstat.geom::volume(L0) * (range(X$df$t)[2] - range(X$df$t)[1])), nX))
    } else {stop(" 'method' argument must be either \"g\" or \"K\ ")}

  }

  lista_0 <- array(NA, dim = c(dim(X_lista)[1:3], k))

  for (i in 1:nX){

    print(i)

    a <- cbind(Z$df, 1)
    b <- cbind(X$df[-i, ], 2)
    colnames(a) <- colnames(b) <- c("x", "y", "t", "id")

    Q_n <- data.frame(rbind(a, b))
    Xi <- X$df[i, ]

    if(verbose) pb <- utils::txtProgressBar(min = 0, max = k, style = 3)
    lista_0 <- sapply(1:k, function(q) {
      if(verbose)  utils::setTxtProgressBar(pb, q)
      if ((inherits(X,"stp")|inherits(Z,"stp"))){
        permutest.stp(perm = perm, Q_n = Q_n, Xi = Xi, method = method)
      } else if ((inherits(X,"stlp")|inherits(Z,"stlp"))){
        permutest.stlp(perm = perm, Q_n = Q_n, Xi = Xi, L0 = L0, method = method)
      }
    }, simplify = "array")
    if(verbose) close(pb)
  }

  mlistas <- apply(lista_0, 1:3, mean, na.rm = TRUE)

  for (i in 1:nX){
    lista_a <- X_lista[, , i]
    lista_b <- mlistas[, , i]
    t2 <- sum((lista_a - lista_b) ^ 2)
    res_test[i] <- t2
    for (j in 1:k){
      lista_a <- lista_0[, , i, j]
      t2_0 <- sum((lista_a - lista_b) ^ 2)
      res_test_d[j,i] <- t2_0
    }
    p[i] <- mean(res_test_d[,i] >= res_test[i])
  }

  id <- which(p <= alpha)

  if ((inherits(X,"stp")|inherits(Z,"stp"))){
    Xsig <- stp(X$df[id, ])
    Xnosig <- stp(X$df[- id, ])
  } else if ((inherits(X,"stlp")|inherits(Z,"stlp"))){
    Xsig <- stp(X$df[id, ], L = L0)
    Xnosig <- stp(X$df[- id, ], L = L0)
  }


  out <- list(p = p, X = X, Z = Z, alpha = alpha, Xsig = Xsig, Xnosig = Xnosig,
              id = id)
  class(out) <- "localtest"
  return(out)
}








