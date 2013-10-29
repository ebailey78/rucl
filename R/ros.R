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
#' @author Eric Bailey

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
  
  #Sort data in decreasing order
  i <- order(x, decreasing = TRUE)
  x <- x[i]
  d <- d[i]
  
  # Detected Values
  xd <- x[d]
  # Empty vector to hold plotting positions for detected values
  pd <- vector("numeric", length(xd))
  # Non-detects
  xnd <- x[!d]
  # Empty vector to hold plotting positions for detected values
  pnd <- vector("numeric", length(xnd))
  
  pmdl <- max(x) + 1     # Previous mdl - starts with dummy value
  pp <- 0                # Previous probability of exceedence - dummy value
  
  # Loop through unique MDLs
  for(i in c(unique(xnd), 0)) {
    
    # Calculate the probability of an exceedence at the current MDL
    x2 <- xd[xd < pmdl]
    mdl.prob <- sum(x2 >= i)/length(x[x< pmdl]) * (1-pp) + pp
    
    # Fill in plotting positions for detected values greater than current MDL
    xi <- xd >= i & xd < pmdl
    pd[xi] <- sapply(seq(sum(xi)), function(j) {
      (1 - pp) - (mdl.prob - pp) / (sum(xi) + 1) * j
    })
    
    # Fill in plotting positions for non-detects equal to current MDL
    xi <- xnd == i
    pnd[xi] <- sapply(seq(sum(xi)), function(j) {
      (1 - mdl.prob) / (sum(xi) + 1) * j
    })
        
    # Remember the current MDL and exceedence probability for next loop
    pp <- mdl.prob
    pmdl <- i
    
  }
  
  val <- rbind(cbind(xd, rep(TRUE, length(xd)), pd), cbind(xnd, rep(FALSE, length(xnd)), pnd))
  val <- val[order(val[, 1], val[, 3], decreasing = TRUE), ]
  
  gp <- gmle(xd) # Estimate gamma parameters
  
  
}