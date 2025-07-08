#' Fit a separable spatio-temporal Poisson process model
#'
#' @param x A \code{stpm} object
#' @param spaceformula A formula for the spatial component. See \link[spatstat.model]{ppm} for details
#' @param timeformula A formula for the temporal component. It fits a log-linear model with the \link{glm} function
#'
#' @return An object of class \code{sepstppm}
#' @export
#'
#' @examples
#' crimesub <- stpm(valenciacrimes$df[101:200, ],
#'            names = colnames(valenciacrimes$df)[-c(1:3)])
#' 
#' mod1 <- sepstppm(crimesub, spaceformula = ~x ,
#'                   timeformula = ~ day)
#'
#'
#'
sepstppm <- function(x, spaceformula, timeformula){

  if (!inherits(x, "stpm")) stop("x should a stpm object")
  
  x0 <- x

  x <- x$df
  
  ot <- x$t

  n <- nrow(x)

  X <- suppressWarnings(spatstat.geom::ppp(x$x, x$y,
                          spatstat.geom::owin(range(x$x), range(x$y))))
  spacemod <- spatstat.model::ppm(as.formula(paste("X", paste(spaceformula, collapse = " "), sep = " ")))
  spaceint <- predict(spacemod, locations = X)
  spaceint_plot <- predict(spacemod)


  stab.expanse1 <- table(x[, all.vars(timeformula)]) |> as.data.frame()
  if(length(all.vars(timeformula))  == 1){
    colnames(stab.expanse1)[1] <- attr(terms(timeformula), "term.labels")
  }
  timemod <- glm(as.formula(paste("Freq", paste(timeformula, collapse = " "), sep = " ")),
                 data = stab.expanse1, family = poisson)

  newdata0 <- data.frame(apply(x[, all.vars(timeformula), drop = F], 2, as.factor))
  timeint <- exp(predict(timemod, newdata = newdata0))
  
  m <- prod(apply(stab.expanse1[ - dim(stab.expanse1)[2]], 2,
                  function(x) length(unique(x))))
  timeint <- n * timeint / m

  stint <- spaceint * timeint #/ n

  out = list(l = stint, x = x0, spaceformula = spaceformula, timeformula = timeformula,
             spacemod = spacemod, timemod = timemod, spaceint = spaceint, timeint = timeint)
  class(out) <- "sepstppm"
  return(out)

}
