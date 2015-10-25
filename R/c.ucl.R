c.ucl <- function(x, d, confidence, N, ...)  {
  
  data.name <- deparse(substitute(x))  # Get the name of the passed variable
  
  # Ensure that there are no NAs in the dataset
  i <- !is.na(x)    
  x <- x[i]
  d <- d[i]
  
  # Various subsets and transformations of the data used by the function
  dx <- x[d]        # Subset of detections only
  dy <- log(dx)     # Natural log of detected values
  n <- length(x)    # Number of samples
  dn <- length(dx)  # Number of samples above detection limit
  
  # These are variables that are needed to calculate other variables and ucl(s) 
  # Perform distribution testing subroutine
  detect.dist <- suppressWarnings(dist.test(dx))  
  # Impute NDs using regression on order (ROS) method (See Helsel or NADA
  # package)
  r <- ros(x, d)                                  
  
  # Calculate gamma maximum liklihood estimate (MLE) of shape, scale, and rate
  # based on detected values
  d.mle <- gamma.mle(dx)              
  # Calcualte gamma maximum liklihood estimate (MLE) of shape, scale, and rate
  # based on gamma ROS imputation
  r.mle <- gamma.mle(r$gamma)        

  # Replace DLs with 1/2 DLs for calculating 1/2 DL UCLs
  dl2n <- x
  dl2n[!d] <- dl2n[!d]/2 
  # Natural log of 1/2 DL dataset
  dl2l <- log(dl2n)
  
  # Dataset censored at maximum detection limit
  dl1 <- x
  dl1[dl1 <= max(x[!d])] = max(x[!d])
  dl1d <- dl1 != max(x[!d])
  
  # Kaplan-Meier estimate of the mean and std error
  km <- ple(x, d)
  # Bootstrap of KM means
  km.boot <- c.bootstrap(x, d, N)
  # Bootstrap of means based on lognormal ROS estimates
  ln.ros.bs <- u.bootstrap(r$ln, N, func = function(x, n, ...) 
    {sum(x)/length(x)})
  
  # Summary Statistics
  ss <- list(n = n, distinct.n = length(unique(dx)), detect.n = dn,
             nondetect.n = n - dn, detect.rate = dn/n, detect.min = min(dx),
             detect.max = max(dx), detect.mean = sum(dx)/dn, 
             detect.sd = sqrt(sum((dx-(sum(dx)/dn))^2)/(dn-1)),
             halfdl.mean = sum(dl2n)/n,
             halfdl.sd = sqrt(sum((dl2n-(sum(dl2n)/n))^2)/(n-1)),
             halfdl.se = sqrt(sum((dl2n-(sum(dl2n)/n))^2)/(n-1)) / sqrt(n),
             norm.ros.mean = sum(r$norm)/n,  
             norm.ros.sd = sqrt(sum((r$norm-(sum(r$norm)/n))^2)/(n-1)), 
             norm.ros.se = sqrt(sum((r$norm-(sum(r$norm)/n))^2)/(n-1)) / sqrt(n), 
             km.mean = km$mean, km.sd = km$se * sqrt(n), km.se = km$se,
             ln.detect.min = min(dy), ln.detect.max = max(dy),
             ln.detect.mean = sum(dy)/dn,
             ln.detect.sd = sqrt(sum((dy-(sum(dy)/dn))^2)/(dn-1)),
             ln.halfdl.mean = sum(dl2l)/n,
             ln.halfdl.sd = sqrt(sum((dl2l-(sum(dl2l)/n))^2)/(n-1)),
             ln.halfdl.se = sqrt(sum((dl2l-(sum(dl2l)/n))^2)/(n-1)) / sqrt(n),
             ln.ros.mean.log = sum(log(r$ln))/n, ln.ros.sd.log = sd(log(r$ln)),
             ln.ros.mean = sum(r$ln)/n,
             ln.ros.sd = sqrt(sum((r$ln-(sum(r$ln)/n))^2)/(n-1)),
             ln.ros.se = sqrt(sum((r$ln-(sum(r$ln)/n))^2)/(n-1))/ sqrt(n),
             g.detect.shape = d.mle$shape, g.detect.scale = d.mle$scale,
             g.detect.rate = d.mle$rate, g.ros.mean = sum(r$gamma)/n,
             g.ros.shape = r.mle$shape, g.ros.scale = r.mle$scale,
             g.ros.rate = r.mle$scale)
  
  # Calculate all UCLs and store in list
  ucls <- list(n.halfdl.t = b.tucl(ss$halfdl.mean, ss$halfdl.se, n, confidence),
               n.ros.t = b.tucl(ss$norm.ros.mean, ss$norm.ros.se, n, confidence),
               n.ros.z = b.zucl(ss$norm.ros.mean, ss$norm.ros.se, confidence),
               ln.halfdl.h = b.hucl(ss$ln.halfdl.mean, ss$ln.halfdl.sd, n, confidence),
               ln.ros.t = b.tucl(ss$ln.ros.mean, ss$ln.ros.se, n, confidence),
               ln.ros.pboot = b.pboot(ln.ros.bs, confidence),
               ln.ros.bcaboot = u.bcaboot(ln.ros.bs, r$ln, ss$ln.ros.mean, n, confidence, N),
               ln.ros.h = b.hucl(ss$ln.ros.mean.log, ss$ln.ros.sd.log, n, confidence),
               g.ros.appgamma = b.appgamma(ss$g.ros.mean, ss$g.ros.shape, n, confidence),
               g.ros.adjgamma = b.adjgamma(ss$g.ros.mean, ss$g.ros.shape, n, confidence),
               o.km.t = b.tucl(km$mean, km$se, n, confidence),
               o.km.z = b.zucl(km$mean, km$se, confidence),
               o.km.tboot = c.tboot(x, d, km$mean, km$se, confidence, N),
               o.km.bcaboot = c.bcaboot(km.boot, x, d, km$mean, n, confidence, N),
               o.km.pboot = b.pboot(km.boot, confidence),
               o.km.cheb95 = b.cheb(km$mean, km$se, 0.95),
               o.km.cheb975 = b.cheb(km$mean, km$se, 0.975),
               o.km.cheb99 = b.cheb(km$mean, km$se, 0.99))
  
  # Calls cu.pick to find recommended UCL(s)
  rec.name <- cu.pick(detect.dist$Recommend, ss$ln.detect.sd, ss$g.detect.shape, n, ss$detect.rate)
  # Pulls the recommended UCL(s) from the ucls list
  rec <- ucls[rec.name]
  
  #Pulls all the pieces into a list that will become a rucl object
  o <- list(type = "censored", confidence = confidence, N = N, data = x, detects = d, 
            data.name = data.name, summary.statistics = ss,  
            dist.test = detect.dist, ros = r, ucls = ucls, rec = rec)
  
  class(o) <- "rucl"
  
  o  
  
}

# A streamlined uncensored ucl function that determines the recommended ucl(s)
# and then only calculates the recommended ucl(s). Returns a numeric vector
# rather than a rucl object and is intended for 'mass production' ucl calcs.
cu.ucl.fast <- function(x, d, confidence, N) {
  
  # Ensure that there are no NAs in the dataset
  i <- !is.na(x)   
  x <- x[i]
  d <- d[i]
  dr <- sum(d)/length(d)  # Detection Rate
  n <- length(x)          # Sample Size
  y <- log(x[d])          # Log transformed detections
  my <- sum(y)/n          # Arithmetic mean of log transformed detections    
  sdy <- sqrt(sum((y-my)^2)/(n-1))  # Standard Deviation of log trans. detects
  
  # Gets the recommended UCL(s)
  rec.name <- cu.pick(suppressWarnings(dist.test(x[d]))$Recommend, 
                      sdy, gshape(x[d]), n, dr)
  
  # Calculate the ple estimate of the mean
  km <- ple(x,d)
  
  # List of unparsed expressions used to calculate recommened UCL(s). These get
  # parsed and evaluated if they are chosen by u.pick().
  ucls <- list(
    o.km.t       = "b.tucl(km$mean, km$se, n, confidence)",
    o.km.pboot   = "b.pboot(c.bootstrap(x, d, N), confidence)",
    o.km.cheb95  = "b.cheb(km$mean, km$se, 0.95)",
    o.km.cheb975 = "b.cheb(km$mean, km$se, 0.975)",
    o.km.cheb99  = "b.cheb(km$mean, km$se, 0.99)",
    o.km.bcaboot = "c.bcaboot(c.bootstrap(x, d, N), x, d, km$mean, n, confidence, N)"
  )
  
  # Loop through the matching UCLs and parse/evaluate the expressions
  v <- sapply(rec.name, function(i) eval(parse(text = ucls[i])))
  names(v) <- rec.name
  v
  
}

# Determine which UCLs may be appropriate for a censored dataset
cu.pick <- function(d, s, k, n, dr) {
  
  # Skewness is based on ln st. dev except for gamma distributions where it is
  # based on the shape parameter
  if(d == "Gamma") s = k
  dr = 1 - dr
  
  # cu.picks is a data.frame representation of the table in the Tech Guide. This
  # dataframe is subsetted by distribution, skewness, and sample size. Then the
  # result is returned as logical vector indicating which UCLs may be valid.
  x <- cu.picks[cu.picks$distribution == d, ]
  x <- x[x$skew.min < s & x$skew.max >= s, ]
  x <- x[x$n.min < n & x$n.max >= n, ]
  x <- x[x$dr.min < dr & x$dr.max >= dr, ]
  x <- x[,8:13]
  if(nrow(x) == 0) {
    o = NULL
  } else {
    o = colnames(x)[x==TRUE]
  }
  
  o
  
}

# A relatively lean and fast bootstrapping function used for censored datasets.
c.bootstrap <- function(x, d, N = as.integer(Sys.getenv("rucl.N")), 
                        func=ple.lite, ...) {
  
  oneRun <- function(N, func, ...) {
    
    x1 <- 1
    d1 <- 1
    
    while(length(unique(x1[d1])) <= 3) {
      i <- sample.int(length(x), replace=TRUE)
      x1 <- x[i]
      d1 <- d[i]
    }
    
    unlist(func(x1, d1, ...))
    
  }
  
  sapply(seq(N), oneRun, func = func, ...=...)
  
} 