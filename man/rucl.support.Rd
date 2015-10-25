\name{Support Functions}
\alias{gscale}
\alias{gshape}
\alias{skew}
\alias{u3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to provide various parameters necessary for UCL calculation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gshape(x, sensitivity = 1e-07, bias.correct = TRUE, neg = "stop")
gscale(x, k, neg = "stop", ...)
skew(x)
u3(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A numeric vector containing the data
  }
  \item{k}{
  An estimate of the the shape parameter
  }
  \item{sensitivity}{
  How precise should the estimate be
  }
  \item{bias.correct}{
  Should bias correction be performed on the result (see details)
  }
  \item{neg}{
  What to do if the MLE calculation encounters a negative value (see details)
  }
  \item{\dots}{
  Additional arguments to be passed to and from other functions
  }
}
\details{
Data that follow a gamma distribution must be positive. When attempting to calculate the maximum likelihood 
estimate of the shape parameter of a dataset containing negative values an error will occur. The 'neg' argument allows the user to choose how to handle this situation. Setting 'neg' to "remove" will remove any value less than or equal to 0 before attempting the calculation. Setting 'neg' to "small" will replace any value less than or equal to 0 with the smallest positve number available on the system. Setting 'neg' to any numeric value will cause values less than or equal to 0 to be replaced with that value. Setting 'neg' to any other value will cause the function to stop with an error (the default behavior).

<There needs to be a brief explanation of the bias correction equation and a reference to the Tech Guide.>


}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eric Bailey <ebailey@idem.in.gov>
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
function (x, k, neg = "stop", ...) 
{
    if (min(x) <= 0) {
        if (neg == "remove") {
            x <- x[x > 0]
        }
        else if (neg == "small") {
            x[x <= 0] <- .Machine$double.eps
        }
        else if (is.numeric(neg)) {
            x[x <= 0] <- neg
        }
        else {
            stop("gamma distributed data cannot contain negative numbers")
        }
    }
    if (missing(k)) 
        k <- gshape(x, ...)
    sum(x)/length(x)/k
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
