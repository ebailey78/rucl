# PLE estimates of variance are much more computationally intensive than for
# the mean. Since most UCLs only require the PLE of the mean this is a faster
# function for use in those situations where the variance isn't needed.
ple.lite <- function(x, d, ...) {
  
  xp <- unique(c(sort(x[d], decreasing = T), min(x)))  
  #xp <- unique(c(sort(x[d], decreasing = T)))  
  
  t <- vector("double", length(xp))
  p <- 1
  
  for(i in 1:(length(t)-1)) {
    
    t1 <- xp[i]
    t2 <- sum(x <= t1)
    p <- p * ((t2 - sum(d[x == t1])) / t2) # (B-C)/B * previous ple
    t[i] <- (t1 - xp[i + 1]) * (1 - p)
    
  }
  
  as.vector(sum(t, xp[length(xp)]))
  
}

# Calculates PLE of mean and variance from a censored dataset
# See 'Improved Methods for Calculating Concentrations Used in Exposure
# Assessments' (BJC/OR-416) for more details.
ple <- function(x, d, ...) {
  
  # Gets list of unique detections plus lowest reading regardless
  d <- as.logical(d)
  xp <- unique(c(min(x), sort(x[d])))  
  
  # Creates a matrix to store intermediate steps of PLE calculation
  t <- matrix(1.0, nrow = length(xp), ncol = 7)
  
  # fills the first row with detected values
  t[, 1] <- xp
  
  # The first row is calculated slightly differently
  t1 <- xp[1]
  t2 <- sum(x <= t1)
  t3 <- sum(d[x == t1])
  t4 <- 0
  t5 <- t1
  t6 <- t5 * t4
  t7 <- t5 - t6
  t[1, ] <- c(t1, t2, t3, t4, t5, t6, t7)
  
  # The last row is also calculated slightly differently
  i <- length(xp)
  t1 <- xp[i]
  t2 <- sum(x <= t1)
  t3 <- sum(d[x == t1])
  t4 <- ((t2 - t3) / t2) 
  t5 <- t1 - t[i - 1, 1]
  t6 <- t5 * t4
  t7 <- t5 - t6
  t[i, ] <- c(t1, t2, t3, t4, t5, t6, t7)
  
  # Loop through rows 2 through n-1 backwards
  for(i in (length(xp) - 1):2) {
    
    t1 <- xp[i]
    t2 <- sum(x <= t1)
    t3 <- sum(d[x == t1])
    t4 <- t[i + 1, 4] * ((t2 - t3) / t2) # (B-C)/B * previous ple
    t5 <- t1 - t[i - 1, 1]
    t6 <- t5 * t4
    t7 <- t5 - t6
    
    t[i, ] <- c(t1, t2, t3, t4, t5, t6, t7)
    
  }
  
  vd <- t[(3 - as.logical(sum(d[x==min(xp)]))):nrow(t), ]
  #  vd <- t[2:nrow(t), ]
  
  # Return a list of the mean and the variance
  list(mean = sum(t[, 7]), 
       se = sqrt(sum(sapply(1:nrow(vd), 
                            function(i, vd) sum(vd[1:i, 6])^2 / (vd[i, 2] * (vd[i, 2] - vd[i, 3])) * vd[i, 3], 
                            vd=vd))
                 * (nrow(t) / (nrow(t) - 1))))
  
}