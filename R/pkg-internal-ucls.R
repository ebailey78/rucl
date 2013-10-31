#### Base functions for calculating UCLs. #### 

#These functions are not intended to be used directly. They are the distilled
#essence of the UCL calculations included in rucl. They assume that all
#necessary inputs have been calculated previously. This approach was taken in an
#attempt to maximize efficiency. Who knows if it really worked though... If you
#are interested in directly accessing the functionality provided by these
#functions please refer to the ucl function.

### Naming Convention ###
# Each of these functions begin with a single character followed by a period.
# This initial character indicated which type of dataset each is intended for as
# follows:
# b - function will work with either a censored or uncensored dataset
# u - function will only work with an uncensored dataset
# c - function will only work with a censored data

### Arguement Convention ### Most UCL calculations require similar inputs. Short
#(usually 1-2 letters) argument names were used for efficiency. The following is a
#key to those short names:
#    m - the arithmetic mean of the dataset (log-transformed for b.hucl)
#   se - the standard error of the mean of the dataset
#   sd - the standard deviation (of the log transformed data)
#    n - the sample size
#  con - the desired confidence for the UCL (1 - alpha) (0.95 by default)
#   u3 - the unbiased third moment estimate
# skew - the skewness of the dataset
#    k - the shape parameter for the gamma distribution
#    v - the variance of the dataset
#    B - a bootstrap dataset containing N bootstrap estimates of a population parameter
#    N - the number of bootstrap iteration that need to be performed
#    x - the dataset being analyzed
#    d - a logical vector indicating whether the values in 'x' are detections or 
#        non-detects: TRUE = detection, FALSE = nondetect

# Student's t UCL
b.tucl <- function(m, se, n, con, ...) m + qt(con, n - 1) * se

# UCL based on central limit theorem
b.zucl <- function(m, se, con, ...) m + qnorm(con) * se

# Modified Student's t UCL
b.modt <- function(m, u3, se, n, con) {
  m + u3 / (6 * (se * sqrt(n))^2 * n) + qt(con, n - 1) * se
}

# Modified Central Limit Theorem (CLT) UCL
b.modz <- function(m, se, skew, n, con) {
  m + (qnorm(con) + (skew / (6 * sqrt(n))) * (1 + 2 * qnorm(con)^2)) * se
}

# Land's H UCL !!PROVIDE WITH LOG TRANSFORMED DATA!!
b.hucl <- function(m, sd, n, con, ...) {
  H <- mns.interpolate(mns.h95, new.x = c(n = n, sd = sd))
  exp(m + 0.5 * sd^2 + sd * H/sqrt(n - 1))
}

# Chebyshev Inequality
b.cheb <- function(m, se, con) m + sqrt((1 / (1 - con)) - 1) * se

# Approximate Gamma
b.appgamma <- function(m, k, n, con, ...) 2 * n * k * m / 
  qchisq(1 - con, 2 * n * k)

# Adjusted Gamma
b.adjgamma <- function(m, k, n, con, ...) {
  b <- max(0, mns.interpolate(mns.adja, new.x = c(n=n, a=1-con)))
  2 * n * k * m / qchisq(b, 2*n*k)
}

# Chebyshev MVUE
b.mvue <- function(m, v, n, con) {
  
  g <- function(t, n) {
    
    w <- 1 
    new.w <- 1
    p <- 1
    
    while(new.w > 5e-10 & p < 100) {
      new.num <- ((n - 1)^(p*2 - 1) * t^p)
      new.den <- (factorial(p) * n^p * prod(seq(n+1, by=2, length.out = p-1)))
      new.w <- new.num/new.den
      w <- w + new.w
      p = p + 1
    }
    
    if(p >= 100) stop("Too many iterations.")
    
    w
    
  }
  
  mean.mvue <- exp(m) * g(v/2, n)
  sd.mvue <- sqrt(exp(2 * m) * (g(v/2, n)^2 - g((n - 2) * v / (n - 1), n)))
  
  mean.mvue + (sqrt((1 / (1 - con)) - 1) * sd.mvue)
  
}

