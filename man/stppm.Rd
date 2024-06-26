% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stppm.R
\name{stppm}
\alias{stppm}
\title{Fit a Poisson process model to a spatio-temporal point pattern}
\usage{
stppm(
  X,
  formula,
  formula_mark = NULL,
  covs = NULL,
  marked = FALSE,
  spatial.cov = FALSE,
  verbose = FALSE,
  mult = 4,
  interp = TRUE,
  parallel = FALSE,
  sites = 1,
  seed = NULL,
  ncube = NULL,
  grid = FALSE,
  ncores = 2,
  lsr = FALSE
)
}
\arguments{
\item{X}{A \code{stp} object}

\item{formula}{An object of class \code{"formula"}:
a symbolic description of the model to be fitted.
The current version only supports formulas depending on the spatial and temporal coordinates:
\code{x}, \code{y}, \code{t}.}

\item{formula_mark}{An object of class \code{"formula"}}

\item{covs}{A list containing \code{stcov} objects of possible spatio-temporal covariates.
It is advisable to construct the \code{stcov} objects with \code{stcov}. 
Each \code{stcov} object should contain the spatio-temporal coordinates and the
covariate values as the fourth column, named as the covariate called in the 
formula.}

\item{marked}{Logical value indicating whether the point process model to be 
fit is multitype. Default to \code{FALSE}.}

\item{spatial.cov}{Logical value indicating whether the point process model to be 
fit depends on spatio-temporal covariates. Default to \code{FALSE}.}

\item{verbose}{Default to \code{FALSE}.}

\item{mult}{The multiplicand of the number of data points,
 for setting the number of dummy
points to generate for the quadrature scheme.}

\item{interp}{Logical value indicating whether to interpolate covariate values 
to dummy points or to use the covariates locations as dummies.
 Default to \code{TRUE}.}

\item{parallel}{Logical values indicating whether to use parallelization to 
interpolate covariates. Default to \code{FALSE}.}

\item{sites}{.....}

\item{seed}{The seed used for the simulation of the dummy points. Default to
\code{NULL}.}

\item{ncube}{Number of cubes used for the cubature scheme.}

\item{grid}{Logical value indicating whether to generate dummy points on a 
regular grid or randomly. Default to \code{FALSE}.}

\item{ncores}{Number of cores to use, if parallelizing. Default to 2.}

\item{lsr}{Logical value indicating whether to use Logistic Spatio-Temporal
Regression or Poisson regression. Default to \code{FALSE}.}
}
\value{
An object of class \code{stppm}. A list of
\describe{
\item{\code{IntCoefs}}{The fitted coefficients}
\item{\code{X}}{The \code{stp} object provided as input}
\item{\code{nX}}{The number of points in \code{X}}
\item{\code{I}}{Vector indicating which points are dummy or data}
\item{\code{y_resp}}{The response variable of the model fitted to the quadrature scheme}
\item{\code{formula}}{The formula provided as input}
\item{\code{l}}{Fitted intensity}
\item{\code{mod_global}}{The \code{glm} object of the model fitted to the quadrature scheme}
\item{\code{newdata}}{The data used to fit the model, without the dummy points}
\item{\code{time}}{Time elapsed to fit the model, in minutes}
}
}
\description{
This function fits a Poisson process model to an observed spatio-temporal
point pattern stored in a \code{stp} object.
}
\details{
We assume that the template model is a Poisson process, with a parametric
intensity or rate function \eqn{\lambda(\textbf{u}, t; \theta)}  with space
and time locations \eqn{\textbf{u} \in
W,  t \in T} and parameters \eqn{\theta \in \Theta.}


Estimation is performed through the fitting of a \code{glm} using a spatio-temporal
version of the quadrature scheme by Berman and Turner (1992).
}
\examples{

set.seed(2)
ph <- rstpp(lambda = 200)
hom1 <- stppm(ph, formula = ~ 1)

## Inhomogeneous
set.seed(2)
pin <- rstpp(lambda = function(x, y, t, a) {exp(a[1] + a[2]*x)}, par = c(2, 6))
inh1 <- stppm(pin, formula = ~ x)

## Inhomogeneous depending on external covariates

set.seed(2)
df1 <- data.frame(runif(100), runif(100), runif(100), rpois(100, 15))
df2 <- data.frame(runif(100), runif(100), runif(100), rpois(100, 15))

obj1 <- stcov(df1, names = "cov1")
obj2 <- stcov(df2, names = "cov2")

covariates <- list(cov1 = obj1, cov2 = obj2)

inh2 <- stppm(pin, formula = ~ x + cov2, covs = covariates, spatial.cov = TRUE)

## Inhomogeneous semiparametric

inh3 <- stppm(pin, formula = ~ s(x, k = 30))

## Multitype

set.seed(2)
dfA <- data.frame(x = runif(100), y = runif(100), t = runif(100), 
                  m1 = rep(c("A"), times = 100))
dfB <- data.frame(x = runif(50), y = runif(50), t = runif(50), 
                  m1 = rep(c("B"), each = 50))

stpm1 <- stpm(rbind(dfA, dfB))

inh4 <- stppm(stpm1, formula = ~ x + s(m1, bs = "re"), marked = TRUE)



}
\references{
Baddeley, A. J., Møller, J., and Waagepetersen, R. (2000). Non-and semi-parametric estimation of interaction in inhomogeneous point patterns. Statistica Neerlandica, 54(3):329–350

Berman, M. and Turner, T. R. (1992). Approximating point process likelihoods with glim. Journal of the Royal Statistical Society: Series C (Applied Statistics), 41(1):31–38

D'Angelo, N., Adelfio, G., and Mateu, J. (2023). Locally weighted minimum contrast estimation for spatio-temporal log-Gaussian Cox processes. Computational Statistics & Data Analysis, 180, 107679.
}
\seealso{
\link{plot.stppm}, \link{print.stppm}, \link{summary.stppm}

\link{locstppm}
}
\author{
Nicoletta D'Angelo and Marco Tarantino
}
