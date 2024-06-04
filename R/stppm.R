#' Fit a Poisson process model to a spatio-temporal point pattern
#'
#' @description
#' This function fits a Poisson process model to an observed spatio-temporal
#' point pattern stored in a \code{stp} object.
#'
#' @details
#' 
#' We assume that the template model is a Poisson process, with a parametric
#' intensity or rate function \eqn{\lambda(\textbf{u}, t; \theta)}  with space
#' and time locations \eqn{\textbf{u} \in
#' W,  t \in T} and parameters \eqn{\theta \in \Theta.}
#'
#'
#' Estimation is performed through the fitting of a \code{glm} using a spatio-temporal
#' version of the quadrature scheme by Berman and Turner (1992).
#'
#' @seealso \link{plot.stppm}, \link{print.stppm}, \link{summary.stppm}
#'
#'
#'
#' @param X A \code{stp} object
#' @param formula An object of class \code{"formula"}:
#' a symbolic description of the model to be fitted.
#' The current version only supports formulas depending on the spatial and temporal coordinates:
#' \code{x}, \code{y}, \code{t}.
#' @param formula_mark An object of class \code{"formula"}
#' @param covs A list containing \code{stcov} objects of possible spatio-temporal covariates.
#' It is advisable to construct the \code{stcov} objects with \code{stcov}. 
#' Each \code{stcov} object should contain the spatio-temporal coordinates and the
#' covariate values as the fourth column, named as the covariate called in the 
#' formula.
#' @param marked Logical value indicating whether the point process model to be 
#' fit is multitype. Default to \code{FALSE}.
#' @param verbose Default to \code{FALSE}.
#' @param spatial.cov Logical value indicating whether the point process model to be 
#' fit depends on spatio-temporal covariates. Default to \code{FALSE}.
#' @param mult The multiplicand of the number of data points,
#'  for setting the number of dummy
#' points to generate for the quadrature scheme.
#' @param interp Logical value indicating whether to interpolate covariate values 
#' to dummy points or to use the covariates locations as dummies.
#'  Default to \code{TRUE}.
#' @param parallel Logical values indicating whether to use parallelization to 
#' interpolate covariates. Default to \code{FALSE}.
#' @param sites .....
#' @param seed The seed used for the simulation of the dummy points. Default to
#'  \code{NULL}.
#' @param ncube Number of cubes used for the cubature scheme.
#' @param grid Logical value indicating whether to generate dummy points on a 
#'  regular grid or randomly. Default to \code{FALSE}.
#' @param ncores Number of cores to use, if parallelizing. Default to 2.
#' @param lsr Logical value indicating whether to use Logistic Spatio-Temporal
#' Regression or Poisson regression. Default to \code{FALSE}.  
#'
#' @return An object of class \code{stppm}. A list of
#' \describe{
#' \item{\code{IntCoefs}}{The fitted coefficients}
#' \item{\code{X}}{The \code{stp} object provided as input}
#' \item{\code{nX}}{The number of points in \code{X}}
#' \item{\code{I}}{Vector indicating which points are dummy or data}
#' \item{\code{y_resp}}{The response variable of the model fitted to the quadrature scheme}
#' \item{\code{formula}}{The formula provided as input}
#' \item{\code{l}}{Fitted intensity}
#' \item{\code{mod_global}}{The \code{glm} object of the model fitted to the quadrature scheme}
#' \item{\code{newdata}}{The data used to fit the model, without the dummy points}
#' \item{\code{time}}{Time elapsed to fit the model, in minutes}
#' }
#'
#'
#' @export
#'
#' @author Nicoletta D'Angelo and Marco Tarantino
#'
#' @seealso
#' \link{locstppm}
#'
#' @examples
#'
#' set.seed(2)
#' ph <- rstpp(lambda = 200)
#' hom1 <- stppm(ph, formula = ~ 1)
#' 
#' ## Inhomogeneous
#' set.seed(2)
#' pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6))
#' inh1 <- stppm(pin, formula = ~ x)
#' 
#' ## Inhomogeneous depending on external covariates
#' 
#' set.seed(2)
#' df1 <- data.frame(runif(100), runif(100), runif(100), rpois(100, 15))
#' df2 <- data.frame(runif(100), runif(100), runif(100), rpois(100, 15))
#'
#' obj1 <- stcov(df1, names = "cov1")
#' obj2 <- stcov(df2, names = "cov2")
#'
#' covariates <- list(cov1 = obj1, cov2 = obj2)
#'
#' inh2 <- stppm(pin, formula = ~ x + cov2, covs = covariates, spatial.cov = TRUE)
#'
#' ## Inhomogeneous semiparametric
#' 
#' inh3 <- stppm(pin, formula = ~ s(x, k = 30))
#' 
#' ## Multitype
#' 
#' set.seed(2)
#' dfA <- data.frame(x = runif(100), y = runif(100), t = runif(100), 
#'                   m1 = rep(c("A"), times = 100))
#' dfB <- data.frame(x = runif(50), y = runif(50), t = runif(50), 
#'                   m1 = rep(c("B"), each = 50))
#' 
#' stpm1 <- stpm(rbind(dfA, dfB))
#' 
#' inh4 <- stppm(stpm1, formula = ~ x + s(m1, bs = "re"), marked = TRUE)
#' 
#' 
#'
#' @references
#' Baddeley, A. J., Møller, J., and Waagepetersen, R. (2000). Non-and semi-parametric estimation of interaction in inhomogeneous point patterns. Statistica Neerlandica, 54(3):329–350
#'
#' Berman, M. and Turner, T. R. (1992). Approximating point process likelihoods with glim. Journal of the Royal Statistical Society: Series C (Applied Statistics), 41(1):31–38
#'
#' D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
#'
stppm <- function(X, formula, formula_mark = NULL, covs = NULL, marked = FALSE, spatial.cov = FALSE,
                  verbose = FALSE, mult = 4, interp = TRUE,
                  parallel = FALSE, sites = 1, seed = NULL, ncube = NULL,
                  grid = FALSE, ncores = 2, lsr = FALSE){
  
  if (!inherits(X, c("stp", "stpm"))) stop("x should be either of class stp or stpm")
  
  time1 <- Sys.time()
  
  if (!is.numeric(mult)) {
    stop("mult should be a numeric value")
  } else {
    if(mult <= 0) {
      stop("mult should be mult > 0")
    }
  } 
  
  if (!is.numeric(sites)) {
    stop("sites should be a numeric value")
  } else {
    if(sites <= 0) {
      stop("sites should be sites > 0")
    }
  } 
  
  if(!is.null(ncube)){
    if(!is.numeric(ncube)) {
    stop("ncube should be a numeric value")
  } else {
    if(ncube <= 0) {
      stop("ncube should be ncube > 0")
    }
  } }
  
  if(is.null(formula_mark)){
    formula_mark <- formula
  }
  
  X0 <- X
  X <- X$df
  
  nX <- nrow(X)
  
  x <- X[,1]
  y <- X[,2]
  t <- X[,3]
  
  if(interp){
    s.region <- splancs::sbox(cbind(x, y), xfrac = 0.01, yfrac = 0.01)
    xr = range(t, na.rm = TRUE)
    xw = diff(xr)
    t.region <- c(xr[1] - 0.01 * xw, xr[2] + 0.01 * xw)
    HomLambda <- nX
    rho <- mult * HomLambda
    
    if (grid == TRUE){
      ff <- floor(rho ^ (1/3))
      x0 <- y0 <- t0 <- 1:ff
      
      scale_to_range <- function(x, new_min, new_max) {
        ((x - min(x)) / (max(x) - min(x))) * (new_max - new_min) + new_min
      }
      
      x0 <- scale_to_range(x0, s.region[1, 1], s.region[2, 1])
      y0 <- scale_to_range(y0, s.region[1, 2], s.region[3, 2])
      t0 <- scale_to_range(t0, t.region[1], t.region[2])
      df0 <- expand.grid(x0, y0, t0)
      dummy_points <- stp(cbind(df0$Var1, df0$Var2, df0$Var3))$df
    } else {
      if(!is.null(seed)) set.seed(seed)
      dummy_points <- rstpp(lambda = rho, nsim = 1, verbose = F,
                            minX = s.region[1, 1], maxX = s.region[2, 1],
                            minY = s.region[1, 2], maxY = s.region[3, 2],
                            minT = t.region[1], maxT = t.region[2])$df
    }
  } else {
    dummy_points <- stp(cbind(covs[[sites]][, 1], covs[[sites]][, 2],
                              covs[[sites]][, 3]))$df
  }
  
  quad_p <- rbind(X[,1:3], dummy_points)
  
  xx <- quad_p[, 1]
  xy <- quad_p[, 2]
  xt <- quad_p[, 3]
  win <- spatstat.geom::box3(xrange = range(xx), yrange = range(xy), zrange = range(xt))
  
  if (is.null(ncube)) {ncube <- .default.ncube(quad_p)}
  ncube <- rep.int(ncube, 3)
  nx <- ncube[1]
  ny <- ncube[2]
  nt <- ncube[3]
  
  nxyt <- nx * ny * nt
  cubevolume <-  spatstat.geom::volume(win) / nxyt  
  volumes <- rep.int(cubevolume, nxyt)
  
  id <- .grid.index(xx, xy, xt, win$xrange, win$yrange, win$zrange, nx, ny, nt)$index
  
  w <- .counting.weights(id, volumes, lsr)
  
  Wdat <- w[1:nX]
  Wdum <- w[-c(1:nX)]
  ndata <- nrow(X)
  ndummy <- nrow(dummy_points)
  
  if(marked == T & spatial.cov == T) {
    
    dati.interpolati <- interp.covariate(X, dummy_points, covs, formula, parallel,
                                         interp, xx, xy, xt, ncores, verbose = verbose)
    
    marked.process <- dummy.marked.result(X, formula = formula_mark, dummy_points, Wdum, Wdat, ndata, ndummy)
    
    z <- c(rep(1, ndata), rep(0, length(marked.process$dumb$data$x)))
    w_final <- c(w[1:ndata], marked.process$Wdumb)
    y_resp <-  if(lsr){z} else {z / w_final}
    
    dati.modello <- data.frame(y_resp, w_final,
                               c(X$x, marked.process$dumb$data$x),
                               c(X$y, marked.process$dumb$data$y),
                               c(X$t, marked.process$dumb$data$z))
    
    colnames(dati.modello) <- c("y_resp", "w", "x", "y", "t")
    
    dati.modello <- cbind(dati.modello, rbind(marked.process$df_marks,
                                              marked.process$total_dummy_marks))
    
    dati.interpolati.rep <- as.data.frame(matrix(NA, nrow = nrow(dati.modello), ncol = ncol(dati.interpolati)-3))
    colnames(dati.interpolati.rep) <- names(dati.interpolati)[-c(1,2,3)]
    q <- 1
    for(i in 4:ncol(dati.interpolati)) {
      dati.interpolati.rep[,q] <- c(dati.interpolati[1:nX,i],
                                    rep(dati.interpolati[(nX+1):nrow(dati.interpolati),i], marked.process$n_comb_levels),
                                    rep(dati.interpolati[1:nX,i], marked.process$n_comb_levels-1))
      q <- q + 1
    }
    dati.cov.marks <- cbind(dati.modello, dati.interpolati.rep)
    
  } else if(marked == F & spatial.cov == T){
    
    dati.interpolati <- interp.covariate(X, dummy_points, covs, formula, parallel,
                                         interp, xx, xy, xt, ncores, verbose = verbose)
    z <- c(rep(1, ndata), rep(0, nrow(dummy_points)))
    y_resp <-  if(lsr){z} else {z / w}
    dati.cov.marks <- cbind(y_resp, w, dati.interpolati)
    
  } else if(marked == F & spatial.cov == F) {
    colnames(quad_p) <- c("x", "y", "t")
    z <- c(rep(1, ndata), rep(0, nrow(dummy_points)))
    y_resp <-  if(lsr){z} else {z / w}
    dati.cov.marks <- cbind(y_resp, w, quad_p)
    
  } else if(marked == T & spatial.cov == F) {
    
    marked.process <- dummy.marked.result(X, formula = formula_mark, dummy_points, Wdum, Wdat, ndata, ndummy)
    
    z <- c(rep(1, ndata), rep(0, length(marked.process$dumb$data$x)))
    w_final <- c(w[1:ndata], marked.process$Wdumb)
    y_resp <-  if(lsr){z} else {z / w_final}
    
    dati.modello <- data.frame(y_resp, w_final,
                               c(X$x, marked.process$dumb$data$x),
                               c(X$y, marked.process$dumb$data$y),
                               c(X$t, marked.process$dumb$data$z))
    
    colnames(dati.modello) <- c("y_resp", "w", "x", "y", "t")
    
    dati.modello <- cbind(dati.modello, rbind(marked.process$df_marks,
                                              marked.process$total_dummy_marks))
    dati.cov.marks <- dati.modello
  }
  
  if(lsr){
    suppressWarnings(mod_global <- try(gam(as.formula(paste("y_resp", paste(formula, collapse = " "),
                                                            "+ offset(log(w))", sep = " ")),
                                           family = binomial,
                                           data = dati.cov.marks), silent = T))
  } else {
    suppressWarnings(mod_global <- try(gam(as.formula(paste("y_resp", paste(formula, collapse = " "), sep = " ")),
                                           weights = w,
                                           family = poisson,
                                           data = dati.cov.marks), silent = T))
  }
  
  
  res_global <- coef(mod_global)
  pred_global <- exp(predict(mod_global, newdata = dati.cov.marks[1:nX, ]))
  
  time2 <- Sys.time()
  
  list.obj <- list(IntCoefs = res_global,
                   X = X0,
                   nX = ndata,
                   I = z,
                   y_resp = y_resp,
                   formula = formula,
                   l = as.vector(pred_global),
                   mod_global = mod_global,
                   newdata = dati.cov.marks[1:ndata, ],
                   ncube = ncube,
                   time = paste0(round(as.numeric(difftime(time1 = time2,
                                                           time2 = time1,
                                                           units = "sec")), 3)," sec"))
  class(list.obj) <- "stppm"
  return(list.obj)
  
}
