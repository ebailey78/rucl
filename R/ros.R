load("C:/R/Packages/rucl - old/data/Oahu.rda")
x <- Oahu[, 1]
d <- Oahu[, 2]

#' Regression on Order Statistics
#' 
#' Impute left-censored datapoints based on normal, lognormal, and gamma 
#' distribution assumptions.
#' 
#' @param x numeric vector containing dataset
#' @param d logical vector indicating detections/non-detects
#' @param \dots additional values to pass to/from other functions
#' @return dataframe containing imputed datasets for normal, lognormal, and
#' gamma assumptions
#' @author Eric Bailey and Nathan Byers

ros <- function(x, d, na.rm = FALSE, ...) {
  
  if(!is.numeric(x)) stop("x must be numeric")
  if(!is.logical(d)) stop("d must be logical")
  if(na.rm) {
    i <- !is.na(x)
    x <- x[i]
    d <- d[i]
  } else if(any(is.na(x))) {
    return(x[FALSE][NA])
  }
    
  p <- vector("numeric", length(x))  # Empty vector to store plotting positions
  pmdl <- max(x) + 1     # Previous mdl - starts with dummy value
  pp <- 0                # Previous probability of exceedence - dummy value
  
  # Loop through unique MDLs
  for(i in c(sort(unique(x[!d]), decreasing = TRUE), 0)) {
    
    # Calculate the probability of an exceedence at the current MDL
    mdl.prob <- sum(x >= i & x < pmdl & d)/length(x[x< pmdl]) * (1-pp) + pp
    
    # Fill in plotting positions for detected values greater than current MDL
    xi <- x >= i & x < pmdl & d
    p[xi] <- sapply(seq(sum(xi)), function(j) {
      (1 - pp) - (mdl.prob - pp) / (sum(xi) + 1) * j
    })
    p[xi][order(x[xi], decreasing = TRUE)] <- p[xi]
    
    # Fill in plotting positions for non-detects equal to current MDL
    xi <- x == i & !d
    p[xi] <- sapply(seq(sum(xi)), function(j) {
      (1 - mdl.prob) / (sum(xi) + 1) * j
    })
        
    # Remember the current MDL and exceedence probability for next loop
    pp <- mdl.prob
    pmdl <- i
    
  }
  
  gp <- gmle(x[d]) # Estimate gamma parameters
  
  qn <- qnorm(p)  # Get standard normal values
  qg <- qgamma(p, shape = gp$shape, scale = gp$scale)  # Standard gamma values
  
  fn <- coef(lm(x[d] ~ qn[d]))  # Fit a linear model (intercept = mean, slope = sd)
  fl <- coef(lm(log(x[d]) ~ qn[d]))  # Fit a linear model (intercept = mean, slope = sd)
  fg <- coef(lm(x[d] ~ qg[d]))  # Fit a linear model (intercept = mean, slope = sd)
  
  pn <- x  # Create predicted values vectors (Normal)
  pl <- x  # (Lognormal)
  pg <- x  # (Gamma)
  pn[!d] <- qn[!d] * fn[2] + fn[1]  # Replace NDs with predictions (Normal)
  pl[!d] <- exp(qn[!d] * fl[2] + fl[1])  # (Lognormal)
  pg[!d] <- qg[!d] * fg[2] + fg[1]       # (Gamma)
  pg[pg <= 0] <- 1e-4
  
  # Create and return a dataframe of the results
  o <- as.data.frame(cbind(x, d, pn, pl, pg))  # Return the data to the original scale
  colnames(o) <- c("original", "detect", "norm", "ln", "gamma")
  o
  
}