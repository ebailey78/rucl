\name{Internal UCLs}
\alias{b.adjgamma}
\alias{b.appgamma}
\alias{b.cheb}
\alias{b.hallboot}
\alias{b.hucl}
\alias{b.modt}
\alias{b.modz}
\alias{b.mvue}
\alias{b.pboot}
\alias{b.tucl}
\alias{b.zboot}
\alias{b.zucl}
\alias{cu.bcaboot}
\alias{cu.tboot}
\alias{uu.bcaboot}
\alias{uu.tboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal rucl UCL functions
}
\description{
Functions used by the 'rucl' package to calculate various UCLs (upper confidence limits of the mean). 
These functions assume that all parameters have been previously calculated. See 'xxxxxxx' for functions intended to be used by the user.
}
\usage{b.adjgamma(m, k, n, a, ...)
b.appgamma(m, k, n, a, ...)
b.cheb(m, se, a)
b.hallboot(x, m, n, skew, se, con, N, ...)
b.hucl(m, sd, n, con, ...)
b.modt(m, u3, se, n, con)
b.modz(m, se, skew, n, con)
b.mvue(m, v, n, a)
b.pboot(B, con)
b.tucl(m, se, n, con, ...)
b.zboot(B, con)
b.zucl(m, se, con, ...)
cu.bcaboot(B, x, d, m, n, con, N, ...)
cu.tboot(x, d, m, se, con, N, ...)
uu.bcaboot(B, x, m, n, con, N, ...)
uu.tboot(x, m, se, con, N, ...)
}

\arguments{
  \item{m}{
the mean of the sample
}
  \item{k}{
an estimate of the gamma distribution shape parameter for the sample
}
  \item{n}{
sample size
}
  \item{a}{
alpha level
}
  \item{se}{
standard error
}
  \item{skew}{
skewness
}
  \item{con}{
confidence level
}
  \item{N}{
number of bootstrap iterations
}
  \item{sd}{
standard deviation
}
  \item{u3}{
%%     ~~Describe \code{a} here~~
}
  \item{v}{
variance
}
  \item{B}{
vector of bootstrap statistics
}
  \item{d}{
dataset used in calculation
}
  \item{\dots}{
additional arguments to be passed to other functions
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a numeric vector of length 1 containing the UCL based on provided parameters.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eric Bailey (Indiana Department of Environmental Management)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m, k, n, a, ...) 
{
    b <- max(0, mns.interpolate(mns.adja, new.x = c(n = n, a = a)))
    2 * n * k * m/qchisq(b, 2 * n * k)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ucls}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
