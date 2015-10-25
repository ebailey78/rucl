\name{rucl Plots}
\alias{plot.rucl}
\alias{barplot.rucl}
\alias{boxplot.rucl}
\alias{hist.rucl}
\alias{qqnorm.rucl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create various plots from rucl objects
}
\description{
Generic functions for producing plots based on rucl objects
}
\usage{
plot.rucl(rucl)
barplot.rucl(rucl)
boxplot.rucl(rucl)
hist.rucl(rucl)
qqnorm.rucl(rucl)
timeplot(rucl, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rucl}{
  A 'rucl' object as created by function ucl
}
}
\details{
These are meant to be called by the generic functions "plot", "boxplot", "barplot", "hist", and "qqnorm".
}
\value{
Calling plot will produce each of the other plots in turn. 
}
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
