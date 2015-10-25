# Performs regression on order estimates based on assumptions of normality,
# lognormality, and gamma distribution.
ros.pp <- function(x, d) {
  
  i <- order(x, decreasing=T) # Sort the data in decreasing order
  x <- x[i]
  d <- d[i]
  
  dl <- c(sort(unique(x[!d]), T),0)  # Get unique DL in decending order
  xp <- vector("double", length(x))  # Create an empty vector to store probabilities
  
  pej <- 0  # Probability of exceedence - Starts at 0 for first bin
  pdl <- max(x) + 1  # Previous Detection Limit - Starts a a value higher than the highest readings
  
  for(j in seq(length(dl))) {
    
    edl <- x[x < pdl & x >= dl[j] & d]  # Values exceeding current detection limit
    A <- length(edl)                    # Same as A in Helsel p. 70
    B <- length(x[x < pdl])             # Same as B in Helsel p. 70
    if(B == 0 & A == 0) {
      C = 0
    } else {
      C = A / B
    }
    n <- length(x[x==dl[j] & !d])       # Number of nondetects at current DL
    ej <- pej + C * (1-pej)         # Current probability estimate
    xp[x %in% edl & d] <- seq(1-pej, 1-ej, length.out=A+2)[-c(1, A+2)]   # Probabilities for values above current DL
    xp[x == dl[j] & !d] <- seq(1-ej, 0, length.out=n+2)[-c(1, n+2)]  # Probabilities for nondetects at current DL
    pdl <- dl[j]  # Set the previous DL for next loop
    pej <- ej     # Set the previous probability estimate for next loop
    
  }
  
  xp
  
}

ros <- function(x, d) {

  i <- order(x, decreasing=T) # Sort the data in decreasing order
  x <- x[i]
  y <- log(x)
  y[y == -Inf] = 2e-16
  d <- d[i]
  
  xp <- ros.pp(x,d)
  
  gmle <- gamma.mle(x[d])

  qn <- qnorm(xp)  # Get the standard normal values based on calculated probabilities
  qg <- qgamma(xp, shape = gmle$shape, scale = gmle$scale)
  
  fn <- coef(lm(x[d] ~ qn[d]))  # Fit a linear model (intercept = mean, slope = sd)
  fl <- coef(lm(y[d] ~ qn[d]))  # Fit a linear model (intercept = mean, slope = sd)
  fg <- coef(lm(x[d] ~ qg[d]))  # Fit a linear model (intercept = mean, slope = sd)
  
  pn <- x  # Create predicted values vectors (Normal)
  pl <- x  # (Lognormal)
  pg <- x  # (Gamma)
  pn[!d] <- qn[!d] * fn[2] + fn[1]  # Replace NDs with predictions (Normal)
  pl[!d] <- exp(qn[!d] * fl[2] + fl[1])  # (Lognormal)
  pg[!d] <- qg[!d] * fg[2] + fg[1]       # (Gamma)
  pg[pg <= 0] <- 1e-4
  #  pg[!d] <- qg[!d]
  
  # Create and return a dataframe of the results
  o <- as.data.frame(cbind(x, d, pn, pl, pg))  # Return the data to the original scale
  colnames(o) <- c("original", "detect", "norm", "ln", "gamma")
  o
  
}