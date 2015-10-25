# Use MLE to estimate shape parameter of gamma distribution
gshape <- function(x, precision = 1e-5, bias.correct = TRUE, neg = "small") {
  
  # Deals with datasets with 0 or negative values
  if(min(x) <= 0) {
    
    if(neg == "remove") {
      x <- x[x > 0]
    } else if(neg == "small") {
      x[x<=0] <- .Machine$double.eps
    } else if (is.numeric(neg)) {
      x[x<=0] <- neg
    } else {
      stop("gamma distributed data cannot contain negative numbers")
    } 
    
  }
  
  n <- length(x)
  m <- sum(x)/n
  y <- log(x)
  
  M <- log(m) - sum(y) / n  # As described on page 47 of Tech Guide
  k = 0                                   # Arbitrary starting point for k
  k1 <- 1/(4*M)*(1 + sqrt(1 + (4/3)*M))   # Starting estimate of k
  c <- 0; stopper <- 100                  # Counter to make sure no inf loops
  
  #Iterative estimatation of k
  while(abs(k-k1) > precision & c < stopper) {
    k <- k1 * ((log(k1) - digamma(k1)) / M)
    k1 = k
    c = c + 1
  }
  
  if(c == stopper) {
    k <- NA
  } else {
    if(bias.correct) k <- (n - 3) * k / n + 2 / (3 * n) # Bias correction (Eq. 2-30)
  }
  
  k
  
}

#Estimate scale parameter of gamma distribution based on shape and mean
gscale <- function(x, k, neg="small", ...) {
  
  if(missing(k)) k <- gshape(x, ...)
  
  sum(x)/length(x)/k
  
}

#Estimates shape, scale and rate of gamma distribution and returns list of values
gamma.mle <- function(x, ...) {
  
  k <- gshape(x, ...)
  t <- gscale(x, k, ...)
  
  list(shape = k, scale = t, rate = 1/t)
  
}