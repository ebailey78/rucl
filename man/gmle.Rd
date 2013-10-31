\name{gmle}
\alias{gmle}
\title{Maximum Likelihood Estimates of gamma Distribution Parameters}
\usage{
  gmle(x, precision = 1e-05, bias.correct = TRUE,
    na.rm = FALSE, neg = "small", ...)
}
\arguments{
  \item{x}{numeric vector containing dataset}

  \item{precision}{indicates to how many decimal places to
  estimate values}

  \item{bias.correct}{\code{logical} indicating whether
  bias correction should be performed}

  \item{neg}{indicates how function should behave if it
  encounters negative values in \code{x}}

  \item{\dots}{additional values to pass to/from other
  functions}
}
\value{
  list containing shape, scale, and rate estimates
}
\description{
  Calculates MLE of shape, scale, and rate parameters for
  the provided dataset assuming a gamma distribution.
}
\examples{
gmle(rgamma(100, shape = 2, scale = 3))
}
\author{
  Eric Bailey
}

