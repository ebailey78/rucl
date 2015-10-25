# The way these tests currently work probably won't pass serious scrutiny. I 
# think they do a reasonable job of detecting the correct distribution as they 
# currently stand but probably don't provide the statistical coverage that they 
# claim to. When a better option for distribution testing presents itself it 
# should be a top priority to replace this module with one that is more 
# statistically valid.
# 
# This package was originally meant to mimic ProUCL as closely as possible, but 
# as development progressed several issues with ProUCL's distribution testing
# arose. Email discussions with the developers resulted in them deciding to
# update the distribution testing in the next version of ProUCL. When this
# version is release the new distribution testing procedures will be evaluated
# and incorporated into this package if possible.

# Shapiro-Wilk Test for normality or lognormality - This is the only test in
# this module that would likely actually pass muster.
sw.test <- function(x, lognormal = FALSE, target.p = 0.05) {
  
  distribution = "Normal"
  # Since this package is mainly intended for environmental applications where
  # negative concentrations would not be possible. It seemed reasonable to
  # replace any 0's or negaitve numbers with a really small number to avoid
  # breaking the log transformation.
  if(lognormal) {
    x[x <= 0] = .Machine$double.eps 
    x = log(x)
    distribution = "Lognormal"
  }
  t <- shapiro.test(x)
  list(method = t$method, 
       distribution = distribution,
       statistic = t$statistic,
       pass = t$p.value >= target.p)
  
}

# Lilliefor's Test for normality or lognormality - This test is likely okay as
# well. !!!!FIND REFERENCE I USED FOR EQUATION ON LINE 57!!!!
lillie.test <- function(x, lognormal = FALSE, target.p = 0.05) {
  
  distribution = "Normal"
  if(lognormal) {
    x[x <= 0] = .Machine$double.eps
    x = log(x)
    distribution = "Lognormal"
  }
  m <- mean(x)
  s <- sd(x)
  t <- ks.test(x, pnorm, mean = m, sd = s)
  n <- length(x)
  D <- t$statistic
  if(n > 100) {
    D <- D * (n/100)^0.49
    n = 100
  }
  
  p <- exp(-7.01256 * D^2 * (n+2.78019) + 2.99587 * D * sqrt(n + 2.78019) - 
             0.122119 + 0.974598/sqrt(n) + 1.67997/n)
  
  list(method = "Lillifors test",
       distribution = distribution,
       statistic = t$statistic,
       pass = as.logical(p >= target.p))
  
}

# Shapiro-Wilkes test for gamma distribution : This is not likely valid. I use 
# pgamma() and qnorm() to transform the dataset to "normal" on the assumption 
# that it was originally gamma with the shape estimated by MLE. I think this is 
# basically valid, the issue comes from estimating the shape because a gamma 
# distributed dataset with a sufficiently high shape parameter (> 20 or so) will
# begin to approximate a normal distribution. Any dataset that is remotely 
# normal-looking will result in an shape mle that will allow the datset to pass 
# for gamma. As such, based on some simulations I performed, it correctly 
# accepts a true gamma-distributed dataset as gamma 95% of the time but accepts 
# non-gamma distrubuted dataset more often than it should. Especially when that 
# dataset is normal. I try to get around this by testing for normality first and
# then gamma. As soon as something better comes along, we should implement it.
ks.gamma <- function(x, target.p = 0.05) {
  
  n <- length(x)
  y <- qnorm(pgamma(x, shape = gshape(x, neg="small")))
  t <- shapiro.test(y)
  if(is.nan(t$p.value)) t$p.value = 0 # Fix for Inf errors with very abhorrent datasets
  list(method = "Shapiro-Wilk Gamma test",
       statistic = t$statistic,
       distribution = "Gamma",
       pass = t$p.value >= target.p)
  
}

# Performs the three tests above and then makes a recommendation on which
# distribution to assume the dataset follows.
dist.test <- function(x) {
  
  x <- x[!is.na(x)]
  if(length(x) > 5000) x <- sample(x, size = 5000)
  
  n.l <- lillie.test(x)$pass
  
  if(min(x) <= 0) {
    l.l <- FALSE
    g.k <- FALSE
  } else {
    l.l <- lillie.test(x, lognormal=T)$pass
    if(gshape(x) > 52) {
      g.k <- FALSE
    } else {
      g.k <- ks.gamma(x)$pass
    }
  }

  if(n.l) {
    dist = "Normal"
  } else if(l.l) {
    dist = "Lognormal"
  } else if(g.k) {
    dist = "Gamma"
  } else {
    dist = "Nonparametric"
  }
  
  list(Normal = n.l, Lognormal = l.l, Gamma = g.k, Recommend = dist)
  
}