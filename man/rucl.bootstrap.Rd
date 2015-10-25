\name{Bootstrap functions}
\alias{cu.bootstrap}
\alias{uu.bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for generating bootstrap samples
}
\description{
Produces bootstrap samples for either censored or uncensored datasets
}
\usage{
cu.bootstrap(x, d, N = as.integer(Sys.getenv("rucl.N")), func = ple.lite, ...)
uu.bootstrap(x, N, func = mean, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector containing the dataset to be bootstrapped
}
  \item{d}{
A logical vector indicating whether readings in x are detections on non-detections, TRUE indicated detection, FALSE indicates non-detection
}
  \item{N}{
The number of bootstrap iterations to be performed
}
  \item{func}{
The function to be bootstrapped, mean by default for uncensored, ple.lite by default for censored datasets
}
  \item{\dots}{
Additional arguments to be passed to and from other functions
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a numeric vector containing the N values calculated from func.}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eric Bailey
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
