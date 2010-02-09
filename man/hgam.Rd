\name{hgam}
\alias{hgam}
\alias{hgam.fit}
\title{
Fitting high-dimensional generalized additive models
}
\description{
\code{hgam} is used to fit high-dimensional generalized additive models.
}
\usage{
hgam(formula, data = NULL, weights, model = LinReg(), 
     nknots = 20, lambda1 = 2, lambda2 = 3, ...)
}
\arguments{
  \item{formula}{
an object of class \code{\link{formula}} (or one that can be coerced to that
class): a symbolic description of the model to be fitted. 
}
  \item{data}{
a data frame.
}
  \item{weights}{ vector of weights. }
  \item{model}{an object of class \code{\link[grplasso]{grpl.model}} implementing the negative
          log-likelihood, gradient, hessian etc. See the documentation
          of \code{\link[grplasso]{grpl.model}} for more details.}
  \item{nknots}{
number of knots.
}
  \item{lambda1}{
grouplasso penalty term.
}
  \item{lambda2}{
smoothing penalty term.
}
  \item{\dots}{
ignored.
}
}
\value{
\code{hgam} returns an object of class \code{hgam}:
  \item{y}{response}
  \item{x}{covariables}
  \item{Btilde}{model matrix}
  \item{coef}{coefficients}
  \item{Btildenew}{function to set up the model matrix for (new) data}
}
\seealso{
\code{\link[grplasso]{grplasso}}
}
\examples{
test.d <- dgp(1000)
test.m <- hgam(y ~ ., data = test.d)
}
\keyword{models}
