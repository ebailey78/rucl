# Function for calculating detailed uncensored ucls
u.ucl <- function(x, confidence, N, ...) {

  data.name <- deparse(substitute(x))  # Get the name of the passed variable
  
  # These are variables that are needed to calculate other variables and ucl(s)
  x <- x[!is.na(x)]               # Ensure that there are no NAs in the dataset
  y <- log(x)                     # Natural log of sample data
  n <- length(x)                  # Number of samples
  m <- sum(x)/n                   # Arithmetic mean of the sample
  sd <- sqrt(sum((x-m)^2)/(n-1))  # Standard deviation of the sample
  se <- sd/sqrt(n)                # Standard error of the mean of the sample
  skew <- skew(x)                 # Skewness of dataset
  
  # Tries to determine the underlying distribution of the dataset
  detect.dist <- suppressWarnings(dist.test(x))
  # Estimates gamma parameters for dataset (shape, scale, rate) using MLE
  gmle <- gamma.mle(x)
  # Bootstrap the arithmetic mean of the dataset
  B <- u.bootstrap(x, N)
  
  # Summary Statistics
  ss <- list(n = n,
             min = min(x), 
             max = max(x), 
             mean = m, 
             geo.mean = exp(sum(y) / n),
             med = median(x), 
             sd = sd, 
             se = se, 
             cov = sd / m,
             skew = skew, 
             ln.min = min(y), 
             ln.max = max(y), 
             ln.mean = sum(y)/n, 
             ln.sd = sd(y),
             g.shape = gmle$shape, 
             g.scale = gmle$scale,
             g.rate = gmle$rate,
             g.mle.mean = gmle$shape * gmle$scale, 
             g.mle.sd = sqrt(gmle$shape * gmle$scale^2)
  )  
  
  # Calculate uncensored UCLs
  ucls <- list(n.tucl = b.tucl(m, se, n, confidence),
               n.modz = b.modz(m, se, skew, n, confidence),
               n.modt = b.modt(m, u3(x), se, n, confidence),
               l.hucl = b.hucl(ss$ln.mean, ss$ln.sd, n, confidence),
               l.mvue95 = b.mvue(ss$ln.mean, ss$ln.sd^2, n, 0.95),
               l.mvue975 = b.mvue(ss$ln.mean, ss$ln.sd^2, n, 0.975),
               l.mvue99 = b.mvue(ss$ln.mean, ss$ln.sd^2, n, 0.99),
               g.appgam = b.appgamma(m, ss$g.shape, n, confidence),
               g.adjgam = b.adjgamma(m, ss$g.shape, n, confidence),
               o.zucl = b.zucl(m, se, confidence),
               o.zboot = b.zboot(B, confidence),
               o.tboot = u.tboot(x, m, se, confidence, N),
               o.hallboot = b.hallboot(x, m, n, ss$skew, se, confidence, N),
               o.pboot = b.pboot(B, confidence),
               o.bcaboot = u.bcaboot(B, x, m, n, confidence, N),
               o.cheb95 = b.cheb(m, se, 0.95),
               o.cheb975 = b.cheb(m, se, 0.975),
               o.cheb99 = b.cheb(m, se, 0.99)
  )
  
  # Calls u.pick to find recommended UCL(s)
  rec.name <- u.pick(detect.dist$Recommend, ss$ln.sd, ss$g.shape, n)
  # Pulls the recommended UCL(s) from the ucls list
  rec <- ucls[rec.name]
  
  #Pulls all the pieces into a list which will become a rucl object
  o <- list(type = "uncensored", confidence = confidence, N = N, 
            data.name = data.name, summary.statistics = ss, 
            dist.test = detect.dist, ucls = ucls, rec = rec)
  
  class(o) <- "rucl"
  
  o
  
}

# A streamlined uncensored ucl function that determines the recommended ucl(s)
# and then only calculates the recommended ucl(s). Returns a numeric vector
# rather than a rucl object and is intended for 'mass production' ucl calcs.
u.ucl.fast <- function(x, confidence, N, ...) {
    
  # Calculates parameters necessary to perform distribution testing
  x <- x[!is.na(x)]
  y <- log(x)
  n <- length(x)
  m <- sum(x) / n
  my <- sum(y) / n
  sd <- sqrt(sum((x - m) ^ 2) / (n - 1))
  sdy <- sqrt(sum((y - my) ^ 2) / (n - 1))
  se <- sd/sqrt(n)

  # Gets the recommended UCL(s)
  rec.name <- u.pick(suppressWarnings(dist.test(x))$Recommend, 
                      sdy, 
                      gshape(x), 
                      n)

  # List of unparsed expressions used to calculate recommened UCL(s). These get
  # parsed and evaluated if they are chosen by u.pick().
  ucls <- list(
    n.tucl     = "b.tucl(m, se, n, confidence)",
    n.modt     = "b.modt(m, u3(x), se, n, confidence)",
    l.hucl     = "b.hucl(my, sdy, n, confidence)",
    o.hallboot = "b.hallboot(x, m, n, skew(x), se, confidence, N)",
    l.mvue95   = "b.mvue(my, sdy^2, n, 0.95)",
    l.mvue975  = "b.mvue(my, sdy^2, n, 0.975)",
    l.mvue99   = "b.mvue(my, sdy^2, n, 0.99)",
    o.cheb95   = "b.cheb(m, se, 0.95)",
    o.cheb975  = "b.cheb(m, se, 0.975)",
    o.cheb99   = "b.cheb(m, se, 0.99)",
    g.appgam   = "b.appgamma(m, gshape(x), n, 1-confidence)",
    o.tboot    = "uu.tboot(x, m, se, confidence, N)",
    g.adjgam   = "b.adjgamma(m, gshape(x), n, 1-confidence)"
  )
  
  v <- sapply(rec.name, function(i) eval(parse(text = ucls[i])))
  names(v) <- rec.name

  v
  
}

# Recommend a UCL based on table in ProUCL tech guide
u.pick <- function(dist, skew, shape, n, ...) {
  
  # Skewness is based on ln st. dev except for gamma distributions where it is
  # based on the shape parameter
  if(dist == "Gamma") skew = shape 
  
  # uu.picks is a data.frame representation of the table in the Tech Guide. This
  # dataframe is subsetted by distribution, skewness, and sample size. Then the
  # result is returned as logical vector indicating which UCLs may be valid.
  x <- uu.picks[uu.picks$distribution == dist, ]
  x <- x[x$skew.min < skew & x$skew.max >= skew, ]
  x <- x[x$n.min < n & x$n.max >= n, ]
  x <- x[, 6:18]
  if(nrow(x) == 0) {
    o = NULL
  } else {
    o = colnames(x)[x==TRUE]
  }
  o  
  
}

# A relatively lean and fast bootstrapping function used for uncensored
# datasets.
u.bootstrap <- function(x, N, func, ...) {
  
  if(missing(func)) func <- function(x) {sum(x)/length(x)}
  
  sapply(seq(N), 
         function(i, x, n, ...) func(sample(x, n, replace=TRUE), ...), 
         x = x, n = length(x), ... = ...)
  
}