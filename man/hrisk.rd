\name{hrisk}
\alias{hrisk}
\title{
Cross-Validation
}
\description{
Cross-validated estimation of the empirical risk for hyper-parameter selection.
}
\usage{
hrisk(object, folds = 10, type = c("cv", "bootstrap", "subsampling"), 
      nlambda1 = 10, lambda2 = 1:10, trace = TRUE, 
      papply = if (require("multicore")) mclapply else lapply)
}
\arguments{
  \item{object}{
an object of class \code{\link{hrisk}}
}
  \item{folds}{
a weight matrix with number of rows equal to the number of observations. The number of columns corresponds to the number of cross-validation runs.
}
  \item{type}{
type of the cross-validation}
  \item{nlambda1}{ignored}
  \item{lambda2}{ignored}
  \item{trace}{ignored}
  \item{papply}{adfa}
}
\details{
If package \code{multicore} is available, \code{hrisk} runs in parallel on cores/processors available.
}
\value{
\code{object} returns an object of class \code{hrisk}.
}
\keyword{models}
