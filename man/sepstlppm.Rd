% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sepstlppm.R
\name{sepstlppm}
\alias{sepstlppm}
\title{Fit a separable spatio-temporal Poisson process model on a linear network}
\usage{
sepstlppm(x, spaceformula, timeformula)
}
\arguments{
\item{x}{A \code{stlpm} object}

\item{spaceformula}{A formula for the spatial component. See \link[spatstat.linnet]{lppm} for details}

\item{timeformula}{A formula for the temporal component. It fits a log-linear model with the \link{glm} function}
}
\value{
An object of class \code{sepstlppm}
}
\description{
Fit a separable spatio-temporal Poisson process model on a linear network
}
\examples{

crimesub <- stpm(valenciacrimes$df[101:200, ],
           names = colnames(valenciacrimes$df)[-c(1:3)],
L = valencianet)

mod1 <- sepstlppm(crimesub, spaceformula = ~x ,
                  timeformula = ~ day)


}
