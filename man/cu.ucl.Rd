\name{cu.ucl}
\alias{cu.ucl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Censored UCLs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cu.ucl(x, d, confidence = as.numeric(Sys.getenv("rucl.confidence")), N = as.numeric(Sys.getenv("rucl.N")))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{confidence}{
%%     ~~Describe \code{confidence} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
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
function (x, d, confidence = as.numeric(Sys.getenv("rucl.confidence")), 
    N = as.numeric(Sys.getenv("rucl.N"))) 
{
    data.name <- deparse(substitute(x))
    i <- !is.na(x)
    x <- x[i]
    d <- d[i]
    n <- length(x)
    dx <- x[d]
    dy <- log(dx)
    dn <- length(dx)
    detect.dist <- suppressWarnings(dist.test(dx))
    r <- ros(x, d)
    d.mle <- gamma.mle(dx)
    r.mle <- gamma.mle(r$gamma)
    dl2n <- x
    dl2n[!d] <- dl2n[!d]/2
    dl2l <- log(dl2n)
    dl1 <- x
    dl1[dl1 <= max(x[!d])] = max(x[!d])
    dl1d <- dl1 != max(x[!d])
    km <- ple(x, d)
    km.boot <- cu.bootstrap(x, d, N)
    ln.ros.bs <- uu.bootstrap(r$ln, N, func = function(x, n, 
        ...) {
        sum(x)/length(x)
    })
    ss.desc <- list("Sample Size", "# of Distinct Detections", 
        "# of Detections", "# of Non-Detects", "Detection Rate", 
        "Minimum Detection", "Maximum Detection", "Arith. Mean of Detections", 
        "SD of Detections", "Arith. Mean using 1/2 DL", "SD using 1/2 DL", 
        "SE using 1/2 DL", "Normal ROS Mean", "Normal ROS SD", 
        "Normal ROS SE", "Kaplan-Meier Mean", "Kaplan-Meier SD", 
        "Kaplan-Meier SE", "Minimum of Log Detections", "Maximum of Log Detections", 
        "Mean of Log Detections", "SD of Log Detections", "Mean of Log Data using 1/2 DL", 
        "SD of Log Data using 1/2 DL", "SE of Log Data using 1/2 DL", 
        "Log ROS-based Mean - Log", "Log ROS-based SD - Log", 
        "Log ROS-based Mean", "Log ROS-based SD", "Log ROS-based SE", 
        "MLE of Gamma Shape - Detects", "MLE of Gamma Scale - Detects", 
        "MLE of Gamma Rate - Detects", "MLE of Gamma Mean - ROS", 
        "MLE of Gamma Shape - ROS", "MLE of Gamma Scale - ROS", 
        "MLE of Gamma Rate - ROS")
    ss <- list(n = n, distinct.n = length(unique(dx)), detect.n = dn, 
        nondetect.n = n - dn, detect.rate = dn/n, detect.min = min(dx), 
        detect.max = max(dx), detect.mean = sum(dx)/dn, detect.sd = sqrt(sum((dx - 
            (sum(dx)/dn))^2)/(dn - 1)), halfdl.mean = sum(dl2n)/n, 
        halfdl.sd = sqrt(sum((dl2n - (sum(dl2n)/n))^2)/(n - 1)), 
        halfdl.se = sqrt(sum((dl2n - (sum(dl2n)/n))^2)/(n - 1))/sqrt(n), 
        norm.ros.mean = sum(r$norm)/n, norm.ros.sd = sqrt(sum((r$norm - 
            (sum(r$norm)/n))^2)/(n - 1)), norm.ros.se = sqrt(sum((r$norm - 
            (sum(r$norm)/n))^2)/(n - 1))/sqrt(n), km.mean = km$mean, 
        km.sd = km$se * sqrt(n), km.se = km$se, ln.detect.min = min(dy), 
        ln.detect.max = max(dy), ln.detect.mean = sum(dy)/dn, 
        ln.detect.sd = sqrt(sum((dy - (sum(dy)/dn))^2)/(dn - 
            1)), ln.halfdl.mean = sum(dl2l)/n, ln.halfdl.sd = sqrt(sum((dl2l - 
            (sum(dl2l)/n))^2)/(n - 1)), ln.halfdl.se = sqrt(sum((dl2l - 
            (sum(dl2l)/n))^2)/(n - 1))/sqrt(n), ln.ros.mean.log = sum(log(r$ln))/n, 
        ln.ros.sd.log = sd(log(r$ln)), ln.ros.mean = sum(r$ln)/n, 
        ln.ros.sd = sqrt(sum((r$ln - (sum(r$ln)/n))^2)/(n - 1)), 
        ln.ros.se = sqrt(sum((r$ln - (sum(r$ln)/n))^2)/(n - 1))/sqrt(n), 
        g.detect.shape = d.mle$shape, g.detect.scale = d.mle$scale, 
        g.detect.rate = d.mle$rate, g.ros.mean = sum(r$gamma)/n, 
        g.ros.shape = r.mle$shape, g.ros.scale = r.mle$scale, 
        g.ros.rate = r.mle$scale)
    ucl.desc <- list("Student's t UCL (1/2 DL)", "Student's t UCL (Normal ROS)", 
        "CLT UCL (Normal ROS)", "Land's H UCL (1/2 DL)", "Student's t UCL (LN ROS)", 
        "Percentile Bootstrap UCL (LN ROS)", "BCA Bootstrap (LN ROS)", 
        "Land's H UCL (LN ROS)", "Approximate Gamma UCL (Gamma ROS)", 
        "Adjusted Gamma UCL (Gamma ROS)", "Student's t UCL (KM)", 
        "CLT UCL (KM)", "Bootstrap-t UCL (KM)", "BCA Percentile Bootstrap UCL (KM)", 
        "Percentile Bootstrap UCL (KM)", "95\% Chebyshev UCL (KM)", 
        "97.5\% Chebyshev UCL (KM)", "99\% Chebyshev UCL (KM)")
    ucls <- list(n.halfdl.t = b.tucl(ss$halfdl.mean, ss$halfdl.se, 
        n, confidence), n.ros.t = b.tucl(ss$norm.ros.mean, ss$norm.ros.se, 
        n, confidence), n.ros.z = b.zucl(ss$norm.ros.mean, ss$norm.ros.se, 
        confidence), ln.halfdl.h = b.hucl(ss$ln.halfdl.mean, 
        ss$ln.halfdl.sd, n, confidence), ln.ros.t = b.tucl(ss$ln.ros.mean, 
        ss$ln.ros.se, n, confidence), ln.ros.pboot = b.pboot(ln.ros.bs, 
        confidence), ln.ros.bcaboot = uu.bcaboot(ln.ros.bs, r$ln, 
        ss$ln.ros.mean, n, confidence, N), ln.ros.h = b.hucl(ss$ln.ros.mean.log, 
        ss$ln.ros.sd.log, n, confidence), g.ros.appgamma = b.appgamma(ss$g.ros.mean, 
        ss$g.ros.shape, n, 1 - confidence), g.ros.adjgamma = b.adjgamma(ss$g.ros.mean, 
        ss$g.ros.shape, n, 1 - confidence), o.km.t = b.tucl(km$mean, 
        km$se, n, confidence), o.km.z = b.zucl(km$mean, km$se, 
        confidence), o.km.tboot = cu.tboot(x, d, km$mean, km$se, 
        confidence, N), o.km.bcaboot = cu.bcaboot(km.boot, x, 
        d, km$mean, n, confidence, N), o.km.pboot = b.pboot(km.boot, 
        confidence), o.km.cheb95 = b.cheb(km$mean, km$se, 0.05), 
        o.km.cheb975 = b.cheb(km$mean, km$se, 0.025), o.km.cheb99 = b.cheb(km$mean, 
            km$se, 0.01))
    rec.name <- uc.pick(detect.dist$Recommend, ss$ln.detect.sd, 
        ss$g.detect.shape, n, ss$detect.rate)
    rec.desc <- ucl.desc[names(ucls) \%in\% rec.name]
    rec <- ucls[rec.name]
    o <- list(type = "censored", confidence = confidence, N = N, 
        data = x, detects = d, data.name = data.name, summary.statistics = ss, 
        ss.desc = ss.desc, distribution.tests = detect.dist, 
        ros = r, ucls = ucls, ucl.desc = ucl.desc, rec = rec, 
        rec.desc = rec.desc)
    class(o) <- "rucl"
    o
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
