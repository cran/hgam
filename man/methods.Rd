\name{methods.hgam}
\alias{methods.hgam}
\alias{print.hgam}
\alias{plot.hgam}
\alias{fitted.hgam}
\alias{coef.hgam}
\alias{predict.hgam}
\alias{logLik.hgam}
\alias{print.hrisk}
\title{
Methods for displaying information about high-dimensional generalized additive models
}
\description{
\code{print.hgam} is used to display some information about the fitted GAMs. 
\code{plot.hgam} plot the marginal effects of the two selected covariates.
\code{fitted.hgam} is used to predict the original data with the fitted GAM model.
\code{coef.hgam} is used to display the fitted coefficients.
\code{predict.hgam} is used to predict (new) data with the fitted GAM model.
\code{logLik.hgam} is used to display the value of the log likelihood function.
\code{print.hrisk} is used to display the results of the hrisk function.
}
\usage{
\method{print}{hgam}(x, \dots)
\method{plot}{hgam}(x, which = NULL, newdata = NULL,
     rug = TRUE, multidim = FALSE, \dots)
\method{fitted}{hgam}(object, \dots)
\method{coef}{hgam}(object, \dots)
\method{predict}{hgam}(object, newdata = NULL, which = NULL, intercept = NULL, \dots)
\method{logLik}{hgam}(object, \dots)
\method{print}{hrisk}(h, \dots)
}
\arguments{
  \item{x}{
an object of class hgam
}
  \item{object}{
an object of class hgam
}
  \item{newdata}{
a data frame. If \code{newdata = NULL} then original data will be used.
}
  \item{which}{
which covariate to be predicted. Character or integer can both be used. If 
\code{which = NULL} all predictors will be used.
}
  \item{multidim}{
when \code{multidim = TRUE} the marginal effects of the two selected covariates is plotted in a three-dimensional grid.
}
  \item{intercept}{
boolean or \code{NULL}. If \code{intercept = NULL} then the function uses the same
intercept options as specified in hgam.
}
  \item{h}{
an object of the class hrisk.
}
  \item{rug}{logical, add rugs to plots.}
  \item{\dots}{
ignored.
}
}
\examples{
test.d <- dgp(1000)
test.m <- hgam(y ~ ., data = test.d)

print(test.m)
coef(test.m)
fitted(test.m)
predict(test.m)
logLik(test.m)
plot(test.m, which = c("x1", "x2"), multidim = TRUE)
}
\keyword{methods}
