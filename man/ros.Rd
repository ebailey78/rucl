\name{ros}
\alias{ros}
\title{Regression on Order Statistics}
\usage{
  ros(x, d, na.rm = FALSE, ...)
}
\arguments{
  \item{x}{numeric vector containing dataset}

  \item{d}{logical vector indicating
  detections/non-detects}

  \item{\dots}{additional values to pass to/from other
  functions}
}
\value{
  dataframe containing imputed datasets for normal,
  lognormal, and gamma assumptions
}
\description{
  Impute left-censored datapoints based on normal,
  lognormal, and gamma distribution assumptions.
}
\author{
  Eric Bailey and Nathan Byers
}

