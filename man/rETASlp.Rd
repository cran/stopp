\name{rETASlp}
\alias{rETASlp}
\title{Simulation of a spatio-temporal ETAS (Epidemic Type Aftershock Sequence) model
on a linear network}
\usage{
rETASlp(
  cat = NULL,
  params = c(0.078915, 0.003696, 0.013362, 1.2, 0.424466, 1.164793),
  betacov = 0.39,
  m0 = 2.5,
  b = 1.0789,
  tmin = 0,
  t.lag = 200,
  xmin = 600,
  xmax = 2200,
  ymin = 4000,
  ymax = 5300,
  iprint = TRUE,
  covdiag = FALSE,
  covsim = FALSE,
  L
)
}
\arguments{
\item{cat}{NULL}

\item{params}{A vector of parameters of the ETAS model to be simulated.
See the 'Details' section.}

\item{betacov}{Numerical array. Parameters of the covariates ETAS model}

\item{m0}{Parameter for the background general intensity of the ETAS model.
In the common seismic analyses it represents the threshold
magnitude.}

\item{b}{1.0789}

\item{tmin}{Minimum value of time.}

\item{t.lag}{200}

\item{xmin}{Minimum of x coordinate range}

\item{xmax}{Maximum of x coordinate range}

\item{ymin}{Minimum of y coordinate range}

\item{ymax}{Maximum of y coordinate range}

\item{iprint}{Default \code{TRUE}}

\item{covdiag}{Default \code{FALSE}}

\item{covsim}{Default \code{FALSE}}

\item{L}{linear network}
}
\value{
A \code{stlp} object
}
\description{
This function simulates a spatio-temporal ETAS
(Epidemic Type Aftershock Sequence) process on a linear network.
It is firstly introduced ans employed for simulation studies in D'Angelo et al. (2021).

It follows the generating scheme for simulating a pattern from an
Epidemic Type Aftershocks-Sequences (ETAS) process
(Ogata and Katsura 1988) with conditional intensity function (CIF) as in
Adelfio and Chiodi (2020), adapted for the space location of events
to be constrained on a linear network.

The simulation on the network is guaranteed by the homogeneous spatial
Poisson processes being generated on the network.

See the 'Details' section.
}
\details{
The CIF of an ETAS
process as in Adelfio and Chiodi (2020) can be written as \deqn{
  \lambda_{\theta}(t,\textbf{u}|\mathcal{H}_t)=\mu f(\textbf{u})+\sum_{t_j<t} \frac{\kappa_0 \exp(\eta_j)}{(t-t_j+c)^p} \{    (\textbf{u}-\textbf{u}_j)^2+d   \}^{-q}  ,
} where

\eqn{\mathcal{H}_t} is the past history of the process up to time
\eqn{t}

\eqn{\mu} is the large-scale general intensity

\eqn{f(\textbf{u})} is
the spatial density

\eqn{\eta_j=\boldsymbol{\beta}' \textbf{Z}_j} is a linear predictor

\eqn{\textbf{Z}_j} the external known covariate vector, including the
magnitude

\eqn{\boldsymbol{\theta}= (\mu, \kappa_0, c, p, d, q, \boldsymbol{\beta})}
are the parameters to be estimated

\eqn{\kappa_0} is a
normalising constant

\eqn{c} and \eqn{p} are characteristic parameters of the
seismic activity of the given region,

 and \eqn{d} and \eqn{q} are two parameters
related to the spatial influence of the mainshock

In the usual ETAS
model for seismic analyses, the only external covariate represents the magnitude,
\eqn{\boldsymbol{\beta}=\alpha}, as
\eqn{\eta_j = \boldsymbol{\beta}' \textbf{Z}_j = \alpha (m_j-m_0)}, where
\eqn{m_j} is the magnitude of the \eqn{j^{th}} event and \eqn{m_0} the threshold
magnitude, that is, the lower bound for which earthquakes with higher
values of magnitude are surely recorded in the catalogue.
}
\examples{

\donttest{
set.seed(5)
X <- rETASlp(cat = NULL, params = c(0.1293688525, 0.003696, 0.013362, 1.2,0.424466,  1.164793),
                         betacov = 0.5, m0 = 2.5, b = 1.0789, tmin = 0, t.lag = 200,
                         xmin = 600, xmax = 2200, ymin = 4000, ymax = 5300, iprint = TRUE,
                         covdiag = FALSE, covsim = FALSE, L = chicagonet)

}

}
\references{
Adelfio, G., and Chiodi, M. (2021). Including covariates in a space-time point process with application to seismicity. Statistical Methods & Applications, 30(3), 947-971.

D???Angelo, N., Adelfio, G., and Mateu, J. (2021). Assessing local differences between the spatio-temporal second-order structure of two point patterns occurring on the same linear network. Spatial Statistics, 45, 100534.

Ogata, Y., and Katsura, K. (1988). Likelihood analysis of spatial inhomogeneity for marked point patterns. Annals of the Institute of Statistical Mathematics, 40(1), 29-39.
}
\author{
Nicoletta D'Angelo and Marcello Chiodi
}
